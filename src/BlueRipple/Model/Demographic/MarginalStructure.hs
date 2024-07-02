{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Strict      #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE UnicodeSyntax #-}
{-# LANGUAGE StandaloneDeriving #-}

module BlueRipple.Model.Demographic.MarginalStructure
  (
    module BlueRipple.Model.Demographic.MarginalStructure
  )
where

import qualified BlueRipple.Model.Demographic.EnrichData as DED

import qualified BlueRipple.Data.Keyed as BRK
import qualified BlueRipple.Data.Types.Demographic as DT


import qualified Control.Foldl as FL
import Control.Lens (Lens', view, over, lens, _2)
import qualified Data.Map.Strict as M
import qualified Data.Map.Merge.Strict as MM
import qualified Data.IntMap.Strict as IM
import qualified Data.Profunctor as PF
import qualified Data.Set as S

import qualified Data.List as List
import qualified Data.Vinyl as V
import qualified Data.Vinyl.TypeLevel as V
import qualified Frames as F

import qualified Flat

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.NLOPT as NLOPT

import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU

import qualified Numeric.ActiveSet as AS

import qualified Data.Vector.Unboxed.Deriving as DU



normalize :: (Functor f, Foldable f) => Lens' a Double -> f a -> f a
normalize l xs = let s = FL.fold (FL.premap (view l) FL.sum) xs in fmap (over l (/ s)) xs

data IsomorphicKeys a b where
  IsomorphicKeys :: (Ord a, BRK.FiniteSet a, Ord b, BRK.FiniteSet b) => (a -> b) -> (b -> a) -> IsomorphicKeys a b

reverseIsomorphism :: IsomorphicKeys a b -> IsomorphicKeys b a
reverseIsomorphism (IsomorphicKeys a2b b2a) = IsomorphicKeys b2a a2b

-- NB: The order of the stencils depends on how we get to k2. Which shouldn't matter as long as
-- we only use them to compute marginals from products to compare to other marginals from
-- products.
reKeyMarginalStructure :: forall k1 k2 w .
                       IsomorphicKeys k2 k1 -> MarginalStructure w k1 -> MarginalStructure w k2
reKeyMarginalStructure ik21 ms = case ik21 of
  (IsomorphicKeys f21 f12) -> case ms of
    MarginalStructure sts ptFld -> MarginalStructure sts' ptFld' where
--      sts' = fmap (expandStencil f21) sts
      sts' = fmap (S.map f12) sts
      ptFld' = PF.dimap (first f21) (sortOn fst . fmap (first f12)) ptFld

-- given a k2 stencil and a map (k1 -> k2), we can generate a k1 stencil by testing each k1 to see if it maps to a k2 which
-- is in the stencil
-- That is, if this is a functor, it's contravariant.
expandStencil :: forall k1 k2 . (BRK.FiniteSet k1, Ord k2, BRK.FiniteSet k2) => (k1 -> k2) -> DED.Stencil Int -> DED.Stencil Int
expandStencil f st = DED.Stencil dil
  where
    rl = S.toList $ BRK.elements @k2
    sis = foldl' (\s i -> flip S.insert s $ rl List.!! i) mempty (DED.stencilIndexes st)
    dl = zip (S.toList $ BRK.elements @k1) [0..]
    dil = snd <$> filter (flip S.member sis . f . fst) dl
{-# INLINEABLE expandStencil #-}


combineMarginalStructures :: forall ok ka kb w .
                             (Ord ok, BRK.FiniteSet ok, Ord ka, BRK.FiniteSet ka, Ord kb, BRK.FiniteSet kb, Monoid w)
                          => Lens' w Double
                          -> (Map ka w -> Map kb w -> Map (ka, kb) w)
                          -> MarginalStructure w (ok, ka) -> MarginalStructure w (ok, kb) -> MarginalStructure w (ok, ka, kb)
combineMarginalStructures wgtLens ip = combineMarginalStructures' wgtLens ip (\(ok, ka, _) -> (ok, ka)) id (\(ok, _, kb) -> (ok, kb)) id id

combineMarginalStructuresF :: forall k ka kb w .
                              (Ord (F.Record k), BRK.FiniteSet (F.Record k)
                              , Ord (F.Record ka), BRK.FiniteSet (F.Record ka)
                              , Ord (F.Record kb), BRK.FiniteSet (F.Record kb)
                              , Ord (F.Record (k V.++ ka V.++ kb))
                              , Ord (F.Record (k V.++ ka)), BRK.FiniteSet (F.Record (k V.++ ka))
                              , Ord (F.Record (k V.++ kb)), BRK.FiniteSet (F.Record (k V.++ kb))
                              , BRK.FiniteSet (F.Record (k V.++ ka V.++ kb))
                              , k F.⊆ (k V.++ ka)
                              , ka F.⊆ (k V.++ ka)
                              , k F.⊆ (k V.++ kb)
                              , kb F.⊆ (k V.++ kb)
                              , (k V.++ ka) F.⊆ (k V.++ ka V.++ kb)
                              , (k V.++ kb) F.⊆ (k V.++ ka V.++ kb)
                              , (k V.++ (ka V.++ kb)) ~ ((k V.++ ka) V.++ kb)
                              , Monoid w
                              )
                           => Lens' w Double
                           -> (Map (F.Record ka) w -> Map (F.Record kb) w -> Map (F.Record ka, F.Record kb) w)
                           -> MarginalStructure w (F.Record (k V.++ ka))
                           -> MarginalStructure w (F.Record (k V.++ kb))
                           -> MarginalStructure w (F.Record (k V.++ ka V.++ kb))
combineMarginalStructuresF wgtLens innerProduct = combineMarginalStructures' wgtLens innerProduct
                                                (F.rcast @(k V.++ ka))
                                                (\r -> (F.rcast @k r, F.rcast @ka r))
                                                (F.rcast @(k V.++ kb))
                                                (\r -> (F.rcast @k r, F.rcast @kb r))
                                          (\(kr, kar, kbr) -> kr F.<+> kar F.<+> kbr)
{-# INLINEABLE combineMarginalStructuresF #-}

combineMarginalStructures' :: forall k ka kb ok ika ikb w .
                             (Ord k, BRK.FiniteSet k
                             , Ord ka, BRK.FiniteSet ka
                             , Ord kb, BRK.FiniteSet kb
                             , Ord ok, BRK.FiniteSet ok
                             , Ord ika, BRK.FiniteSet ika
                             , Ord ikb, BRK.FiniteSet ikb
                             , Monoid w)
                           => Lens' w Double
                           -> (Map ika w -> Map ikb w -> Map (ika, ikb) w)
                           -> (k -> ka)
                           -> (ka -> (ok, ika))
                           -> (k -> kb)
                           -> (kb -> (ok, ikb))
                           -> ((ok, ika, ikb) -> k)
                           -> MarginalStructure w ka
                           -> MarginalStructure w kb
                           -> MarginalStructure w k
combineMarginalStructures' wgtLens outerProduct kaF expandA kbF expandB collapse msA msB = MarginalStructure subsets prodFld
  where
--    expandedStencilsA = fmap (expandStencil kaF) $ msStencils msA
--    expandedStencilsB = fmap (expandStencil kbF) $ msStencils msB
--    sts = expandedStencilsA <> expandedStencilsB
    subsetsA = fmap (mapSubset kaF) $ msSubsets msA
    subsetsB = fmap (mapSubset kbF) $ msSubsets msB
    subsets = subsetsA <> subsetsB

    aMapTblFld :: FL.Fold (k, w) (Map ok (Map ika w))
    aMapTblFld = FL.fold (normalizedTableMapFld wgtLens) <$> (fmap (fmap (first expandA)) $ FL.premap (first kaF) $ msProdFld msA)
    bMapTblFld :: FL.Fold (k, w) (Map ok (Map ikb w))
    bMapTblFld = FL.fold (normalizedTableMapFld wgtLens) <$> (fmap (fmap (first expandB)) $ FL.premap (first kbF) $ msProdFld msB)
    prodFld :: FL.Fold (k, w) [(k, w)]
    prodFld = fmap (fmap (first collapse)) $ tableProductL outerProduct <$> aMapTblFld <*> bMapTblFld
{-# INLINEABLE combineMarginalStructures' #-}

mapSubset :: (Ord k, BRK.FiniteSet k, Ord k', BRK.FiniteSet k') => (k' -> k) -> Set k -> Set k'
mapSubset f = S.fromList . mconcat . fmap inverseImage . S.toList
  where
    inverseImage k = filter ((== k) . f) $ S.toList BRK.elements

identityMarginalStructure :: forall k w . (Ord k, BRK.FiniteSet k, Monoid w)
                          => Lens' w Double -> MarginalStructure w k
identityMarginalStructure wgtLens = MarginalStructure eachAsSubset (M.toList <$> normalizeAndFillMapFld wgtLens)
  where
    eachAsSubset = fmap S.singleton $ S.toList BRK.elements
{-# INLINEABLE identityMarginalStructure #-}


-- NB: The stencil order is unspecified. So this only works if the stencils are used as a map from a joint distribution to a marginal one
data MarginalStructure w k where
  MarginalStructure :: (Monoid w, BRK.FiniteSet k, Ord k) => [Set k] -> FL.Fold (k, w) [(k, w)]-> MarginalStructure w k

subsetToStencil :: (BRK.FiniteSet k, Ord k) => Set k -> DED.Stencil Int
subsetToStencil s = DED.Stencil $ fmap snd $ filter fst $ zip (M.elems comboMap) [0..]
  where
    allMap = M.fromList $ fmap (, False) $ S.toList BRK.elements
    subsetMap = M.fromList $ fmap (, True) $ S.toList s
    comboMap = M.unionWith (||) allMap subsetMap

stencilToSubset :: (BRK.FiniteSet k, Ord k) => DED.Stencil Int -> Set k
stencilToSubset (DED.Stencil indices) = S.fromList $ fmap snd $ filter fst $ zip (IM.elems comboMap) (S.toList elts)
  where
    elts = BRK.elements
    nElts = S.size elts
    allMap = IM.fromList $ fmap (, False) [0..(nElts - 1)]
    subsetMap = IM.fromList $ fmap (, True) $ indices
    comboMap = IM.unionWith (||) allMap subsetMap

msSubsets :: MarginalStructure w k -> [Set k]
msSubsets (MarginalStructure subsets _) = subsets
{-# INLINE msSubsets #-}

msStencils :: MarginalStructure w k -> [DED.Stencil Int]
msStencils (MarginalStructure subsets _) = fmap subsetToStencil subsets
{-# INLINE msStencils #-}

msProdFld :: MarginalStructure w k -> FL.Fold (k, w) [(k, w)]
msProdFld (MarginalStructure _ ptF) = ptF
{-# INLINE msProdFld #-}

msNumCategories :: forall k w . MarginalStructure w k -> Int
msNumCategories ms = case ms of
  MarginalStructure _ _ -> S.size $ BRK.elements @k
{-# INLINEABLE msNumCategories #-}


marginalProductFromJointFld :: Lens' w Double -> MarginalStructure w k -> FL.Fold (k, w) [(k, w)]
marginalProductFromJointFld wgtLens ms =
  let nFld = FL.premap (view $ _2 . wgtLens) FL.sum
      f n normalizedPT = fmap (over (_2 . wgtLens) (* n)) normalizedPT
  in case ms of
    MarginalStructure _ ptFld -> f <$> nFld <*> ptFld

constMap :: (BRK.FiniteSet k, Ord k) => a -> Map k a
constMap x = M.fromList $ fmap (, x) $ S.toList BRK.elements
{-# INLINEABLE constMap #-}

zeroMap :: (Monoid w, BRK.FiniteSet k, Ord k) => Map k w
zeroMap = constMap mempty
{-# INLINEABLE zeroMap #-}
{-# SPECIALIZE zeroMap :: (BRK.FiniteSet k, Ord k, Num b) => Map k (Sum b) #-}
{-# SPECIALIZE zeroMap :: (BRK.FiniteSet k, Ord k) => Map k CellWithDensity #-}

summedMap :: (Monoid b, Ord a) => FL.Fold (a, b) (Map a b)
summedMap = FL.foldByKeyMap FL.mconcat --fmap getSum <$> FL.premap (second Sum) appendingMap
{-# INLINE summedMap #-}
{-# SPECIALIZE summedMap :: (Ord a, Num b) => FL.Fold (a, Sum b) (Map a (Sum b)) #-}
{-# SPECIALIZE summedMap :: Ord a => FL.Fold (a, CellWithDensity ) (Map a CellWithDensity) #-}

zeroFillSummedMapFld :: (BRK.FiniteSet k, Ord k, Monoid w) => FL.Fold (k, w) (Map k w)
zeroFillSummedMapFld = fmap (<> zeroMap) summedMap
{-# INLINEABLE zeroFillSummedMapFld #-}
{-# SPECIALIZE zeroFillSummedMapFld ::  (BRK.FiniteSet k, Ord k, Num b) => FL.Fold (k, Sum b) (Map k (Sum b)) #-}
{-# SPECIALIZE zeroFillSummedMapFld ::  (BRK.FiniteSet k, Ord k) => FL.Fold (k, CellWithDensity) (Map k CellWithDensity) #-}

--normalized :: (Foldable f, Functor f, Normalizable c) => f Double -> f Double
--normalized xs = let s = FL.fold FL.sum xs in fmap ( / s) xs

normalizeMapFld :: (Ord k, Monoid w) =>  Lens' w Double -> FL.Fold (k, w) (Map k w)
normalizeMapFld wgtLens = normalize wgtLens <$> summedMap
{-# INLINEABLE normalizeMapFld #-}

normalizeAndFillMapFld :: (BRK.FiniteSet k, Ord k, Monoid w)
                       => Lens' w Double -> FL.Fold (k, w) (Map k w)
normalizeAndFillMapFld wgtLens = normalize wgtLens <$> zeroFillSummedMapFld
{-# INLINEABLE normalizeAndFillMapFld #-}

unNestTableProduct :: Map outerK (Map (a, b) x) -> [((outerK, a, b), x)]
unNestTableProduct = concatMap (\(ok, mab) -> fmap (\((a, b), x) -> ((ok, a, b), x)) $ M.toList mab) . M.toList
{-# INLINEABLE unNestTableProduct #-}

-- Nested maps make the products easier
-- this fills in missing items with 0
tableMapFld :: forall outerK x w . (BRK.FiniteSet outerK, Ord outerK, BRK.FiniteSet x, Ord x, Monoid w) => FL.Fold ((outerK, x), w) (Map outerK (Map x w))
tableMapFld = FL.premap keyPreMap $ fmap (<> constMap zeroMap) (FL.foldByKeyMap zeroFillSummedMapFld) where
    keyPreMap ((o, k), n) = (o, (k, n))
{-# INLINEABLE tableMapFld #-}


--mapWgtLens :: (w -> Double, (Double -> Double) -> w -> w) -> (Map x w -> Double, (Double -> Double) -> Map x w -> Map x w)
--mapWgtLens (getWgt, modWgt) = (FL.fold (FL.premap getWgt FL.sum), fmap . modWgt)
{-
mapWgtLens' :: Ord x => Traversal w w Double Double -> Traversal (Map x w) (Map x w) Double Double
mapWgtLens' l afa mxw = s2
  where
    s1 = fmap (\w -> (w, afa $ view l w)) mxw
    s2 = traverse (\(w, fx) -> fmap (\x -> set l x w) fx) s1
-}
normalizedTableMapFld :: forall outerK x w . (BRK.FiniteSet outerK, Ord outerK, BRK.FiniteSet x, Ord x, Monoid w)
                      =>  Lens' w Double -> FL.Fold ((outerK, x), w) (Map outerK (Map x w))
normalizedTableMapFld wgtLens = FL.premap keyPreMap $ fmap (normalizeMap . (<> constMap zeroMap)) (FL.foldByKeyMap zeroFillSummedMapFld)
  where
    keyPreMap ((o, k), n) = (o, (k, n))
    sum' mm = FL.fold FL.sum $ fmap (FL.fold (FL.premap (view wgtLens) FL.sum)) mm
    normalizeMap mm = let mmSum = sum' mm in fmap (fmap (over wgtLens (/ mmSum))) mm
{-# INLINEABLE normalizedTableMapFld #-}


tableProductL' ::  forall outerK a b w .
                   (Ord outerK, BRK.FiniteSet a, Ord a, BRK.FiniteSet b, Ord b, Monoid w)
               => (Map a w -> Map b w -> Map (a, b) w)
               -> Map outerK (Map a w)
               -> Map outerK (Map b w)
               -> [w]
tableProductL' innerProduct ma mb = snd <$> tableProductL innerProduct ma mb
{-# INLINEABLE tableProductL' #-}

tableProductL ::  forall outerK a b w .
                  (Ord outerK, BRK.FiniteSet a, Ord a, BRK.FiniteSet b, Ord b, Monoid w)
              => (Map a w -> Map b w -> Map (a, b) w)
              -> Map outerK (Map a w)
              -> Map outerK (Map b w)
              -> [((outerK, a, b), w)]
tableProductL innerProduct ma mb = unNestTableProduct $ tableProduct innerProduct ma mb
{-# INLINEABLE tableProductL #-}

-- This is not symmetric. Table a is used for weights on outer key
-- Also, this treats missing entries as 0
tableProduct :: forall outerK a b w .
                (Ord outerK, BRK.FiniteSet a, Ord a, BRK.FiniteSet b, Ord b, Monoid w)
             => (Map a w -> Map b w -> Map (a, b) w)
             -> Map outerK (Map a w)
             -> Map outerK (Map b w)
             -> Map outerK (Map (a, b) w)
tableProduct outerProduct aTableMap bTableMap = MM.merge
                                                (MM.mapMissing whenMissing)
                                                (MM.mapMissing whenMissing)
                                                (MM.zipWithMatched whenMatched)
                                                aTableMap bTableMap
  where
    whenMissing _ _ = zeroMap
    whenMatched _ = outerProduct
{-# INLINEABLE tableProduct #-}



innerProductSum :: (Ord a, Ord b, Eq x, Fractional x) => Map a (Sum x) -> Map b (Sum x) -> Map (a, b) (Sum x)
innerProductSum ma mb = M.fromList [f ae be | ae <- M.toList ma, be <- fracs mb]
      where
        f (a, x) (b, y) = ((a, b), Sum $ getSum x * getSum y)
        fracs :: (Eq x, Fractional x) => Map z (Sum x) -> [(z, Sum x)]
        fracs m = let s = getSum (FL.fold FL.mconcat m) in M.toList $ if s == 0 then m else fmap (Sum . (/ s) . getSum) m
{-# INLINEABLE innerProductSum #-}

data CellWithDensity = CellWithDensity { cwdWgt :: !Double, cwdDensity :: !Double } deriving (Show, Generic, Flat.Flat)

updateWgt :: CellWithDensity -> Double -> CellWithDensity
updateWgt (CellWithDensity _ wd) x = CellWithDensity x wd
{-# INLINE updateWgt #-}

cwdWgtLens :: Lens' CellWithDensity Double --(CellWithDensity -> Double, (Double -> Double) -> CellWithDensity -> CellWithDensity)
cwdWgtLens = lens cwdWgt updateWgt
{-# INLINEABLE cwdWgtLens #-}

safeDiv :: Double -> Double -> Double
safeDiv x y = if y /= 0 then x / y else 0

instance Semigroup CellWithDensity where
  (CellWithDensity x xwd) <> (CellWithDensity y ywd) =
    CellWithDensity
    (x + y)
    ((x * xwd + y * ywd) `safeDiv` (x + y))
  {-# INLINEABLE (<>) #-}

instance Semigroup CellWithDensity => Monoid CellWithDensity where
  mempty = CellWithDensity 0 0
  mappend = (<>)

{-
updateWeightCWD :: Double -> CellWithDensity -> CellWithDensity
updateWeightCWD x (CellWithDensity wgt pwDensity) = CellWithDensity x (pwDensity * if wgt == 0 then 1 else x / wgt)
-}


data DensityProduct = GMDensity
                    | SimpleDensity
                    | SVDDensity
                    | PosConstrainedDensity deriving (Show, Eq)

configuredInnerProductCWD :: (Ord a, BRK.FiniteSet a, Ord b, BRK.FiniteSet b
                             , Show a, Show b)
                          => DensityProduct ->  Map a CellWithDensity -> Map b CellWithDensity -> Map (a, b) CellWithDensity
configuredInnerProductCWD densP ma mb = f ma mb where
  f =  case densP of
         GMDensity -> innerProductGM
         SimpleDensity -> innerProductCWD
         SVDDensity -> innerProductCWD'
         PosConstrainedDensity -> innerProductCWD''

cwdToRec :: CellWithDensity -> F.Record [DT.PopCount, DT.PWPopPerSqMile]
cwdToRec (CellWithDensity p d) = round p F.&: d F.&: V.RNil

innerProductGM :: (Ord a, Ord b) => Map a CellWithDensity -> Map b CellWithDensity -> Map (a, b) CellWithDensity
innerProductGM ma mb = M.fromList [f ae be | ae <- M.toList ma, be <- forProd mb]
  where
    f (a, CellWithDensity xw _) (b, (yw', _)) = ((a, b), CellWithDensity (xw * yw') gmd)
    gmdF = fmap snd $ DT.densityAndPopFld' DT.Geometric (const 1) cwdWgt cwdDensity
    gmdA = FL.fold gmdF ma
    gmdB = FL.fold gmdF mb
    gmd = sqrt $ gmdA * gmdB
    forProd :: Map x CellWithDensity -> [(x, (Double, Double))]
    forProd m =
      let evenFrac = 1 / realToFrac (M.size m)
          sumWgtFld = FL.premap cwdWgt FL.sum
--          sumWgtdDensityFld = FL.premap (\t -> cwdWgt t * cwdDensity t) FL.sum
          sumWgts = FL.fold sumWgtFld m
          g (CellWithDensity x xwd)
            | sumWgts == 0 = (evenFrac, xwd)
            | otherwise = (x / sumWgts, xwd)
      in M.toList $ fmap g m
{-# INLINEABLE innerProductGM #-}

innerProductCWD :: (Ord a, Ord b) => Map a CellWithDensity -> Map b CellWithDensity -> Map (a, b) CellWithDensity
innerProductCWD ma mb = M.fromList [f ae be | ae <- M.toList ma, be <- forProd mb]
  where
    f (a, CellWithDensity xw xwd) (b, (yw', ywd')) = ((a, b), CellWithDensity (xw * yw') (xwd * ywd'))
    forProd :: Map x CellWithDensity -> [(x, (Double, Double))]
    forProd m =
      let evenFrac = 1 / realToFrac (M.size m)
          sumWgtFld = FL.premap cwdWgt FL.sum
          sumWgtdDensityFld = FL.premap (\t -> cwdWgt t * cwdDensity t) FL.sum
          (sumWgts, sumWgtdDensities) = FL.fold ((,) <$> sumWgtFld <*> sumWgtdDensityFld) m
          g (CellWithDensity x xwd)
            | sumWgts == 0 = (evenFrac, evenFrac)
            | sumWgtdDensities == 0 = (x / sumWgts, evenFrac)
            | otherwise = (x / sumWgts, x * xwd / sumWgtdDensities)
      in M.toList $ fmap g m
{-# INLINEABLE innerProductCWD #-}

innerProductCWD' :: (BRK.FiniteSet a, Ord a, BRK.FiniteSet b, Ord b, Show a, Show b) => Map a CellWithDensity -> Map b CellWithDensity -> Map (a, b) CellWithDensity
innerProductCWD' ma mb =
  let la = FL.fold FL.length ma
      lb = FL.fold FL.length mb
      na = let x = FL.fold (FL.premap cwdWgt FL.sum) ma in x -- trace (show ma) x
      nb = let x = FL.fold (FL.premap cwdWgt FL.sum) mb in x -- trace (show mb) x -- this should be the same and maybe we shoudl check?
      ma' = let f = if na == 0 then (const $ 1 / realToFrac la) else (/ na) in fmap (over cwdWgtLens f) ma
      mb' = let f = if nb == 0 then (const $ 1 / realToFrac lb) else (/ nb) in fmap (over cwdWgtLens f) mb
      countAB ((a, wa), (b, wb)) = ((a, b), realToFrac na * cwdWgt wa * cwdWgt wb)
      gmDens ((_, wa), (_, wb)) = sqrt $ cwdDensity wa * cwdDensity wb
      ab = [(a, b) | a <- M.toList ma', b <- M.toList mb']
      countsAB = fmap countAB ab
      gmDensV = VS.fromList $ fmap gmDens ab
      vecA x =  VS.fromList $ fmap (\((a, _), (_, bw)) -> if x == a then cwdWgt bw else 0) ab
      vecB x = VS.fromList $ fmap (\((_, aw), (b, _)) -> if x == b then cwdWgt aw else 0) ab
      rows = (fmap vecA $ M.keys ma') <> (fmap vecB $ M.keys mb')
      m = let x = LA.fromRows rows in x --trace ("m=" ++ toString (LA.dispf 2 x)) x
      rhsV = let x = (VS.fromList $ fmap cwdDensity (M.elems ma') <> fmap cwdDensity (M.elems mb')) - m LA.#> gmDensV in x
--             in trace ("rhs=" ++ toString (DED.prettyVector x)) x
      solM = LA.linearSolveSVD m $ LA.fromColumns [rhsV]
      solV = let x = List.head $ LA.toColumns $ solM in x -- trace ("sol" ++ toString (DED.prettyVector x)) x
      solL = VS.toList $ solV + gmDensV
      mkAssoc (k, wgt) d = (k, CellWithDensity wgt d)
  in M.fromList $ zipWith mkAssoc countsAB solL
{-# INLINEABLE innerProductCWD' #-}

innerProductCWD'' :: (BRK.FiniteSet a, Ord a, BRK.FiniteSet b, Ord b, Show a, Show b) => Map a CellWithDensity -> Map b CellWithDensity -> Map (a, b) CellWithDensity
innerProductCWD'' ma mb =
  let --la = FL.fold FL.length ma
      --lb = FL.fold FL.length mb
      na = let x = FL.fold (FL.premap cwdWgt FL.sum) ma in x -- trace (show ma) x
      nb = let x = FL.fold (FL.premap cwdWgt FL.sum) mb in x -- trace (show mb) x -- this should be the same and maybe we shoudl check?
      ma' = if na == 0 then ma else fmap (over cwdWgtLens (/ na)) ma
      mb' = if nb == 0 then mb else fmap (over cwdWgtLens (/ nb)) mb
      probAB ((a, p), (b, q)) = ((a, b), cwdWgt p * cwdWgt q)
      ab = [(a, b) | a <- M.toList ma', b <- M.toList mb']
      probsAB = fmap probAB ab
      vecA x =  VS.fromList $ fmap (\((a, _), (_, q)) -> if x == a then cwdWgt q else 0) ab
      vecB x = VS.fromList $ fmap (\((_, p), (b, _)) -> if x == b then cwdWgt p else 0) ab
      rows = (fmap vecA $ M.keys ma') <> (fmap vecB $ M.keys mb')
      m = let x = LA.fromRows rows in x --trace ("m=" ++ toString (LA.dispf 2 x)) x
      rhsV = let x = VS.fromList $ fmap cwdDensity (M.elems ma') <> fmap cwdDensity (M.elems mb') in x
      gmDens ((_, wa), (_, wb)) = sqrt $ cwdDensity wa * cwdDensity wb
      gmDensV = VS.fromList $ fmap gmDens ab
      solL = positiveConstrainedSolve m rhsV gmDensV
{-
--             in trace ("rhs=" ++ toString (DED.prettyVector x)) x
      solM = LA.linearSolveSVD m $ LA.fromColumns [rhsV]
      solV = let x = List.head $ LA.toColumns $ solM in x -- trace ("sol" ++ toString (DED.prettyVector x)) x
      solL = VS.toList solV
-}
      mkAssoc n (k, p) d = (k, CellWithDensity (n * p) d)
      allZero = M.fromList $ fmap (\(k, _) -> (k, CellWithDensity 0 0)) probsAB
  in if na == 0 || nb == 0 then allZero else M.fromList $ zipWith (mkAssoc na) probsAB $ VS.toList solL
{-# INLINEABLE innerProductCWD'' #-}

{-
nnlsSolve :: LA.Matrix Double -> LA.Vector Double -> LA.Vector Double -> K.Sem r (LA.Vector Double)
nnlsSolve m rhs = do
  let asConfig = AS.defaultConfig {AS.cfgStart = StartGS }
  AS.optimalNNLS (K.logLE K.Info) asConfig m rhs
-}

positiveConstrainedSolve :: LA.Matrix Double -> LA.Vector Double -> LA.Vector Double -> LA.Vector Double
positiveConstrainedSolve m rhs v0 =  if fst (objective v0) < 1 then v0 else g where
  mTm = LA.tr m LA.<> m
  vLen = VS.length v0
  objective v = (LA.norm_2 $ (m LA.#> v) - rhs, 2 * mTm LA.#> v)
--      constraintD n v = (negate $ v VS.! n, (LA.assoc vLen 0 [(n, negate 1)] :: LA.Vector Double))
--      nlConstraintsD cTol = fmap (\n -> NLOPT.InequalityConstraint (NLOPT.Scalar (constraintD n)) cTol) [0..(vLen - 1)]
  maxIters = 1000
  relTol = 1
--  relTolV = VS.replicate vLen relTol
  nlStop = NLOPT.ParameterRelativeTolerance relTol :| [NLOPT.MaximumEvaluations maxIters, NLOPT.MinimumValue 1]
  nlAlgo = NLOPT.SLSQP objective [NLOPT.LowerBounds $ VS.replicate vLen 0] [] []
  nlProblem = NLOPT.LocalProblem (fromIntegral vLen) nlStop nlAlgo
  nlSol = NLOPT.minimizeLocal nlProblem v0
  g = case nlSol of
    Left _result -> let x = v0 in trace ("positiveConstrainedSolve failed. Returning v0") x
    Right solution -> case NLOPT.solutionResult solution of
      NLOPT.MAXEVAL_REACHED -> let x = NLOPT.solutionParams solution in trace ("positiveConstrainedSolve: NLOPT Solver hit max evaluations (" <> show maxIters
                                                                               <> "). rhs=" <> toString (DED.prettyVector rhs) <> "; v="
                                                                                <> toString (DED.prettyVector x)
                                                                                <> "; m=" <> toString (LA.dispf 2 m)
                                                                                <> "; mv - rhs=" <> toString (DED.prettyVector $ (m LA.#> x) - rhs)
                                                                                <> "; ||mv - rhs||=" <> show (LA.norm_2 $ (m LA.#> x) - rhs)
                                                                              )
                                                                        x
      NLOPT.MAXTIME_REACHED -> let x = NLOPT.solutionParams solution in trace ("positiveConstrainedSolve: NLOPT Solver hit max time.") x
      _ -> NLOPT.solutionParams solution


DU.derivingUnbox "CellWithDensity"
  [t|CellWithDensity -> (Double, Double)|]
  [|\(CellWithDensity w d) -> (w, d)|]
  [|\(w, d) -> CellWithDensity w d|]

{-
-- given a k2 product fold and a map (k1 -> k2) we can generate a k1 product zero-info product fold by assuming every k1 which maps to a k2
-- gets an equal share of that k2s weight
-- This will always fold to Just something because k1 is a finite set. But I can't get rid of the Maybe.
expandProdFld :: forall k1 k2 . (Ord k2, BRK.FiniteSet k1)
              => (k1 -> k2) -> FL.Fold (k2, Double) [(k2, Double)] -> FL.Fold (k1, Double) [(k1, Double)]
expandProdFld f k2Fld = expandList <$> fldA
  where
    fldA = FL.premap (first f) k2Fld
    k2k1Map = foldl' (\m k1 -> M.insertWith (<>) (f k1) [k1] m) mempty BRK.elements
    eF (k2, x) = case M.lookup k2 k2k1Map of
      Just k1s -> fmap (, x / realToFrac (length k1s)) k1s
      Nothing -> []
    expandList :: [(k2, Double)] -> [(k1, Double)]
    expandList = concatMap eF
{-# INLINEABLE expandProdFld #-}


expandMarginalStructure :: (Ord k1, BRK.FiniteSet k1) => (k1 -> k2) -> MarginalStructure k2 -> MarginalStructure k1
expandMarginalStructure f m2 = case m2 of
  MarginalStructure sts2 ptFld2 -> MarginalStructure (fmap (expandStencil f) sts2) (expandProdFld f ptFld2)
{-# INLINEABLE expandMarginalStructure #-}
-}
