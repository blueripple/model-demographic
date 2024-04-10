{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE UnicodeSyntax #-}
{-# LANGUAGE StandaloneDeriving #-}

module BlueRipple.Model.Demographic.NullSpaceBasis
  (
    module BlueRipple.Model.Demographic.NullSpaceBasis
  )
where

import qualified BlueRipple.Model.Demographic.EnrichData as DED
import qualified Data.Map.Strict as M
import qualified Data.List as List
import qualified Data.Set as Set
import qualified Data.Text as T

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Array as A
import qualified Numeric.LinearAlgebra.Array.Util as A

import qualified Data.Vector.Storable as VS

-- # of dimensions at each index
newtype Dimensions = Dimensions { dimensions :: [Int] }

-- # a set of indices into the array
newtype Indices = Indices { indices :: [Int] }

-- a subset specified by choosing some indices
newtype Subset = Subset { subset :: [Int] } deriving newtype (Eq, Ord)

data Known a where
  KnownMarginal :: a -> Known a -- we know the powerset of this
  KnownSubset :: a -> Known a -- we know just this
   deriving (Functor, Foldable, Traversable)
   deriving stock (Show)

knownSubsets :: Known Subset -> Set.Set Subset
knownSubsets (KnownMarginal s) = Set.fromList $ fmap Subset $ powerset $ subset s
knownSubsets (KnownSubset s) = Set.singleton s

knownsToSubsets :: Traversable f => f (Known Subset) -> Set.Set Subset
knownsToSubsets = Set.unions . fmap knownSubsets

knownCharText :: Known [Char] -> T.Text
knownCharText (KnownMarginal cs) = "!" <> T.pack cs
knownCharText (KnownSubset cs) = "+" <> T.pack cs

-- given a set of stencils, i.e., rows of constraints on distributions,
-- which are possibly linearly-dependent,
-- compute an orthonormal basis for them
-- as well as an orthonormal basis for their orthogonal complement, the null space
nullSpacePartitionSVD :: Int -> [DED.Stencil Int] -> (LA.Matrix LA.R, LA.Matrix LA.R)
nullSpacePartitionSVD n sts =
  let cM = DED.mMatrix n sts
      (_, s, v') = LA.svd cM
      k = LA.ranksv (2.2e-16) (max (LA.rows cM) (LA.cols cM)) (LA.toList s)
      nVs = LA.tr $ LA.dropColumns k v'
      nonRedundantConstraints = LA.tr $ LA.takeColumns k v'
  in (nonRedundantConstraints, nVs)

contrastBasis :: Int -> LA.Matrix LA.R
contrastBasis n = LA.fromRows (rowOnes : negIRows)
  where
    rowOnes = LA.fromList $ List.replicate (n-1) 1
    negIRows = LA.toRows $ negate $ LA.ident (n-1)

contrastBasis' :: Int -> LA.Matrix LA.R
contrastBasis' n = LA.fromRows (iRows <> [rowOnes])
  where
    rowOnes = LA.fromList $ List.replicate (n-1) (-1)
    iRows = LA.toRows $ LA.ident (n-1)


-- matrix of rows of 1s or appropriate slices of the contrast basis
aVecs :: Dimensions -> Subset -> Indices -> [LA.Vector LA.R]
aVecs dims sub subsetIs =
  let sWi = M.fromList $ zip (subset sub) (indices subsetIs)
      row (dimPos, dimSize) = case M.lookup dimPos sWi of
        Nothing -> LA.fromList $ List.replicate dimSize 1
        Just sElt -> List.head $ LA.toRows $ LA.dropRows (sElt - 1) (LA.tr $ contrastBasis' dimSize)
  in fmap row $ zip [1..] (dimensions dims)

infixl 9 #
(#) :: [Int] -> [Double] -> A.Array Double
ds # cs = A.listArray ds cs :: A.Array Double

rArray :: Char -> LA.Vector LA.R -> A.Array LA.R
rArray i m =
  let n = LA.size m
  in [n] # (LA.toList m) A.! [i]

asStackedRowVec :: A.Array LA.R -> LA.Vector LA.R
asStackedRowVec x = List.head $ LA.toColumns $ A.matrixator x (reverse $ fmap A.iName $ A.dims x) []

asStackedColVec :: VS.Storable LA.R => A.Array LA.R -> LA.Vector LA.R
asStackedColVec x = List.head $ LA.toColumns $ A.matrixator x (fmap A.iName $ A.dims x) []

-- given as a subset of bs, compute the vector corresponding to
interactionV :: Dimensions -> Subset -> Indices -> LA.Vector LA.R
interactionV dims sub subsetIs =
  let a = aVecs dims sub subsetIs
      nDims = length (dimensions dims)
      indexLabels = take nDims ['a'..]
      rT k = rArray (indexLabels List.!! (k - 1)) (a List.!! (k - 1))
      rTs = fmap rT [2..nDims]
  in asStackedRowVec $ List.foldl' (A.|*|) (rT 1) rTs

interactionVC :: Dimensions -> Subset -> Indices -> LA.Vector LA.R
interactionVC dims sub subsetIs =
  let a = aVecs dims sub subsetIs
      nDims = length (dimensions dims)
      indexLabels = take nDims ['a'..]
      rT k = rArray (indexLabels List.!! (k - 1)) (a List.!! (k - 1))
      rTs = fmap rT [2..nDims]
  in asStackedColVec $ List.foldl' (A.|*|) (rT 1) rTs


subsetIndices :: Dimensions -> Subset -> [Indices]
subsetIndices dims subs = fmap Indices $ traverse (\si -> [1..((dimensions dims List.!! (si - 1)) - 1)]) $ subset subs

subsetInteractionBasis :: Dimensions -> Subset -> LA.Matrix LA.R
subsetInteractionBasis dims sub = LA.fromColumns $ fmap oneVec subIs
  where
    oneVec = interactionV dims sub
    subIs = subsetIndices dims sub

{-
powerSetInteractionBasis :: Dimensions -> Subset -> LA.Matrix LA.R
powerSetInteractionBasis dims subs =
  LA.fromColumns
  $ mconcat
  $ fmap (LA.toColumns . subsetInteractionBasis dims) subsets
  where
    subsets = fmap Subset $ powerset $ subset subs

interactionBasis :: Dimensions -> Set.Set Subset -> LA.Matrix LA.R
interactionBasis dims subs =
  LA.fromColumns
  $ mconcat
  $ fmap (LA.toColumns . subsetInteractionBasis dims) subsets
  where
    uniques = Set.toList . Set.fromList
    subsets = fmap Subset
              $ uniques
              $ mconcat
              $ fmap (powerset . subset) $ Set.toList subs

-}
-- k is a phantom here to at least assure we are using a category map with the correct types
-- But this is upsettingly Stringy
data CatMap k where
  CatMap :: [Char] -> (Char -> Maybe Int) -> (Char -> Maybe Int) -> CatMap k

catMapDimensions :: CatMap k -> Maybe Dimensions
catMapDimensions (CatMap cns lp ls) = Dimensions <$> traverse ls cns

reTypeCatMap :: forall k' k . CatMap k -> CatMap k'
reTypeCatMap (CatMap c ls lp) = CatMap c ls lp

data CatsAndKnowns k = CatsAndKnowns { camCatMap :: CatMap k, camKnowns :: [Known [Char]]}

reTypeCatsAndKnowns :: forall k' k . CatsAndKnowns k -> CatsAndKnowns k'
reTypeCatsAndKnowns (CatsAndKnowns cm ms) = CatsAndKnowns (reTypeCatMap cm) ms

knownsText :: CatsAndKnowns k -> T.Text
knownsText cam = T.intercalate "_" $ fmap knownCharText $ camKnowns cam

catsText :: CatsAndKnowns k -> T.Text
catsText cam = let (CatMap c _ _) = camCatMap cam in T.pack c

mkCatMap :: forall k . [(Char, Int)] -> CatMap k
mkCatMap catNamesAndSizes = CatMap allCats lookupPos lookupSize where
  allCats = fst <$> catNamesAndSizes
  lookupPos c = M.lookup c $ M.fromList $ zip (fst <$> catNamesAndSizes) [1..]
  lookupSize c = M.lookup c $ M.fromList catNamesAndSizes

catSubset :: CatMap k -> [Char] -> Maybe [Int]
catSubset (CatMap _ f _) = traverse f

powerset :: [a] -> [[a]]
powerset l = [] : go l
  where
    go [] = []
    go (x : xs) = let ps = go xs in [x] : (fmap (x :) ps) <> ps

subsetFromKnownChars :: CatMap k -> Known [Char] -> Maybe (Known Subset)
subsetFromKnownChars (CatMap _ lp _) = traverse (fmap Subset . traverse lp)

interactionBasisCM' :: Traversable f => CatMap k -> f (Known [Char]) -> Maybe (LA.Matrix LA.R)
interactionBasisCM' cm knownSubsetsC = do
  dims <- catMapDimensions cm
  knownSubsets <- knownsToSubsets <$> traverse (subsetFromKnownChars cm) knownSubsetsC
  pure $ LA.fromColumns $ mconcat $ fmap (LA.toColumns . subsetInteractionBasis dims) $ Set.toList knownSubsets
--  knownSubsets <- Set.fromList <$> (traverse subsetFromChars knownSubsets)
--  pure $ interactionBasis dims subsets

interactionBasisCM :: CatsAndKnowns k -> Maybe (LA.Matrix LA.R)
interactionBasisCM cam = interactionBasisCM' (camCatMap cam) (camKnowns cam)

nullSpacePartitionCM' :: Traversable f => CatMap k -> f (Known [Char]) -> Maybe (LA.Matrix LA.R, LA.Matrix LA.R)
nullSpacePartitionCM' cm knownSubsetsC = do
  dims <- catMapDimensions cm
  knownSubsets <- knownsToSubsets <$> traverse (subsetFromKnownChars cm) knownSubsetsC
  let allSubsets = knownsToSubsets [KnownMarginal $ Subset [1..(length $ dimensions dims)]]
--  powerSetOfKnownSubsets <- Set.fromList . fmap Subset . mconcat . fmap (powerset . subset) <$> (traverse subsetFromChars $ Set.toList knownSubsets)
--  unionOfKnownSubsets <- Set.fromList <$> (traverse subsetFromChars $ Set.toList knownSubsets)
  constraints <- interactionBasisCM' cm knownSubsetsC
  let nullSpace = LA.fromColumns $ mconcat $ fmap (LA.toColumns . subsetInteractionBasis dims) $ Set.toList $ Set.difference allSubsets knownSubsets
  pure (LA.tr constraints, LA.tr nullSpace)

nullSpacePartitionCM :: CatsAndKnowns k -> Maybe (LA.Matrix LA.R, LA.Matrix LA.R)
nullSpacePartitionCM cam = nullSpacePartitionCM' (camCatMap cam) (camKnowns cam)
