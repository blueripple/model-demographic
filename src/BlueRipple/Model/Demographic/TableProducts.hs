{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
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

module BlueRipple.Model.Demographic.TableProducts
  (
    module BlueRipple.Model.Demographic.TableProducts
  , module Numeric.ActiveSet
  )
where

import qualified BlueRipple.Model.Demographic.EnrichData as DED
import qualified BlueRipple.Model.Demographic.MarginalStructure as DMS
import qualified BlueRipple.Model.Demographic.NullSpaceBasis as DNS

import qualified BlueRipple.Data.Keyed as BRK

import qualified Knit.Report as K
import qualified Polysemy.Error as PE

import qualified Control.MapReduce.Simple as MR
import qualified Frames.MapReduce as FMR
import qualified Frames.Streamly.InCore as FSI

import qualified Numeric.ActiveSet as AS
import Numeric.ActiveSet (ActiveSetConfiguration(..), defaultActiveSetConfig, EqualityConstrainedSolver(..), Logging(..), NNLS_Start(..), precomputeFromE)
import qualified Control.Foldl as FL
import qualified Data.List as L
import qualified Data.Map.Strict as M
import qualified Data.Map.Merge.Strict as MM
import qualified Data.Set as S

import qualified Data.List as List
import qualified Data.Vinyl as V
import qualified Data.Vinyl.TypeLevel as V
import qualified Frames as F
import qualified Frames.Melt as F
import qualified Frames.Transform as FT
import qualified Numeric as Numeric
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.NLOPT as NLOPT
import qualified Data.Vector.Storable as VS

import Control.Lens (Lens', view, lens)

import qualified Flat
import Flat.Instances.Vector()

normalizedVec :: VS.Vector Double -> VS.Vector Double
normalizedVec v = VS.map (/ VS.sum v) v

stencils ::  forall a b.
             (Ord a
             , Ord b
             , BRK.FiniteSet b
             )
         => (b -> a) -> [DED.Stencil Int]
stencils bFromA = M.elems $ FL.fold (DED.subgroupStencils bFromA) $ S.toList BRK.elements
{-# INLINEABLE stencils #-}


zeroKnowledgeTable :: forall outerK a . (Ord outerK, BRK.FiniteSet outerK, Ord a, BRK.FiniteSet a) => Map outerK (Map a Double)
zeroKnowledgeTable = DMS.constMap (DMS.constMap (1 / num)) where
    num = realToFrac $ S.size (BRK.elements @outerK) * S.size (BRK.elements @a)

-- does a pure product
-- for each outerK we split each bs cell into as many cells as are represented by as
-- and we split in proportion to the counts in as
frameTableProduct :: forall outerK as bs count r .
                     (V.KnownField count
                     , DED.EnrichDataEffects r
                     , (bs V.++ (outerK V.++ as V.++ '[count])) ~ (((bs V.++ outerK) V.++ as) V.++ '[count])
                     , outerK F.⊆ (outerK V.++ as V.++ '[count])
                     , outerK F.⊆ (outerK V.++ bs V.++ '[count])
                     , bs F.⊆ (outerK V.++ bs V.++ '[count])
                     , FSI.RecVec (outerK V.++ as V.++ '[count])
                     , FSI.RecVec (bs V.++ (outerK V.++ as V.++ '[count]))
                     , F.ElemOf (outerK V.++ as V.++ '[count]) count
                     , F.ElemOf (outerK V.++ bs V.++ '[count]) count
                     , Show (F.Record outerK)
                     , BRK.FiniteSet (F.Record bs)
                     , Ord (F.Record bs)
                     , Ord (F.Record outerK)
                     , V.Snd count ~ Int
                     )
                  => F.FrameRec (outerK V.++ as V.++ '[count])
                  -> F.FrameRec (outerK V.++ bs V.++ '[count])
                  -> K.Sem r (F.FrameRec (bs V.++ outerK V.++ as V.++ '[count]))
frameTableProduct base splitUsing = DED.enrichFrameFromModel @count (fmap (DED.mapSplitModel round realToFrac) . splitFLookup) base
  where
    splitFLookup = FL.fold (DED.splitFLookupFld (F.rcast @outerK) (F.rcast @bs) (realToFrac @Int @Double . F.rgetField @count)) splitUsing

-- to get a set of counts from the product counts and null-space weights
-- take the weights and left-multiply them with the vectors
-- to get the weighted sum of those vectors and then add it to the
-- product counts


-- k is a phantom here, unmentioned in the data. But it forces us to line things up with the marginal structure, etc.
-- Given N entries in table and a marginal structure with a Q dimensional null space
-- Given C constraints, first matrix is the C x N matrix of constraints coming from the structure of known marginals
-- Second matrix, P, is projections, Q x N, each row is a basis vector of the null space.
-- So Pv = projection of v onto the null space
-- For PCAW version:
-- vector is covariance of each Whitened projections
-- Third matrix, R, is rotation from null space basis to eigenvector basis. Columns are eigenvectors.
-- So R'P projects onto the null space and rotates into the basis of eigenvectors of covariance matrix
-- And P'R rotates a vector in the nullspace expressed in the basis of eigenvectors of the covariance matrix and
-- rotates back to the SVD computed null space and then injects into the space of the full table.

data NullVectorProjections k where
  NullVectorProjections :: (Ord k, BRK.FiniteSet k) =>  DNS.NullSpacePartition -> NullVectorProjections k
  PCAWNullVectorProjections :: (Ord k, BRK.FiniteSet k) =>  DNS.NullSpacePartition -> LA.Vector Double -> LA.Matrix Double -> NullVectorProjections k

nullSpacePartition :: NullVectorProjections k -> DNS.NullSpacePartition
nullSpacePartition (NullVectorProjections x) = x
nullSpacePartition (PCAWNullVectorProjections x _ _) = x

nvpConstraints :: NullVectorProjections k -> LA.Matrix Double
nvpConstraints  = DNS.knownM . nullSpacePartition

nvpProj :: NullVectorProjections k -> LA.Matrix Double
nvpProj = DNS.unknownM . nullSpacePartition

nvpRot :: NullVectorProjections k -> LA.Matrix Double
nvpRot (NullVectorProjections nsp) = LA.ident $ LA.rows $ DNS.unknownM nsp
nvpRot (PCAWNullVectorProjections _ _ r) = r

nvpToEither :: NullVectorProjections k
             -> Either
             DNS.NullSpacePartition
             (DNS.NullSpacePartition, LA.Vector Double, [LA.Vector Double])
nvpToEither (NullVectorProjections nsp) = Left nsp
nvpToEither (PCAWNullVectorProjections nsp s r) = Right (nsp, s, LA.toRows r)

nvpFromEither :: (Ord k, BRK.FiniteSet k)
              => Either
                 DNS.NullSpacePartition
                 (DNS.NullSpacePartition, LA.Vector Double, [LA.Vector Double])
              -> NullVectorProjections k
nvpFromEither e = case e of
  Left nsp -> NullVectorProjections nsp
  Right (nsp, s, rr) -> PCAWNullVectorProjections nsp s (LA.fromRows rr)

instance (Ord k, BRK.FiniteSet k) => Flat.Flat (NullVectorProjections k) where
  size  = Flat.size . nvpToEither
  encode = Flat.encode . nvpToEither
  decode = nvpFromEither <$> Flat.decode


{-
instance (Ord k, BRK.FiniteSet k) => Flat.Flat (NullVectorProjections k) where
  size (NullVectorProjections c p r) = Flat.size (LA.toRows c, LA.toRows p, LA.toRows r)
  encode (NullVectorProjections c p r) = Flat.encode (LA.toRows c, LA.toRows p, LA.toRows r)
  decode = (\(cVs, pVs, rVs) -> NullVectorProjections (LA.fromRows cVs) (LA.fromRows pVs) (LA.fromRows rVs)) <$> Flat.decode
-}


data ProjectionsToDiff k where
  RawDiff :: NullVectorProjections k -> ProjectionsToDiff k
  AvgdDiff :: LA.Matrix Double -> NullVectorProjections k -> ProjectionsToDiff k

type PTDEither k = Either (NullVectorProjections k) ([LA.Vector Double], NullVectorProjections k)

instance Flat.Flat (NullVectorProjections k) => Flat.Flat (ProjectionsToDiff k) where
  size (RawDiff nvp) = Flat.size @(PTDEither k) (Left nvp)
  size (AvgdDiff m nvp) = Flat.size @(PTDEither k) (Right (LA.toRows m, nvp))
  encode (RawDiff nvp) = Flat.encode @(PTDEither k) (Left nvp)
  encode (AvgdDiff m nvp) = Flat.encode @(PTDEither k) (Right (LA.toRows m, nvp))
  decode = f <$> Flat.decode where
    f = \case
      Left nvp -> RawDiff nvp
      Right (mRows, nvp) -> AvgdDiff (LA.fromRows mRows) nvp

nullVectorProjections :: ProjectionsToDiff k -> NullVectorProjections k
nullVectorProjections (RawDiff x) = x
nullVectorProjections (AvgdDiff _ x) = x

numProjections :: NullVectorProjections k -> Int
numProjections = fst . LA.size .  DNS.unknownM . nullSpacePartition

fullToProjM :: NullVectorProjections k -> LA.Matrix LA.R
fullToProjM (NullVectorProjections nsp) = DNS.unknownM nsp
fullToProjM (PCAWNullVectorProjections nsp _ rM) = LA.tr rM LA.<> DNS.unknownM nsp

projToFullM :: NullVectorProjections k -> LA.Matrix LA.R
projToFullM (NullVectorProjections nsp) = DNS.ipUnknownM nsp
projToFullM (PCAWNullVectorProjections nsp _ rM) = DNS.ipUnknownM nsp <> rM

projToFull :: NullVectorProjections k -> LA.Vector LA.R -> LA.Vector LA.R
projToFull nvps v = projToFullM nvps LA.#> v

fullToProj :: NullVectorProjections k -> LA.Vector LA.R -> LA.Vector LA.R
fullToProj nvps v = fullToProjM nvps LA.#> v

projToNullM :: NullVectorProjections k -> LA.Matrix LA.R
projToNullM nvps = projToFullM nvps LA.<> fullToProjM  nvps

projToNonNullM :: NullVectorProjections k -> LA.Matrix LA.R
projToNonNullM nvps = let f2p = fullToProjM nvps in LA.ident (LA.rows f2p) - (projToFullM nvps LA.<> f2p)


{-
baseNullVectorProjections :: forall w k . (BRK.FiniteSet k) => DMS.MarginalStructure w k -> NullVectorProjections k
baseNullVectorProjections ms = case ms of
  DMS.MarginalStructure _ _ -> NullVectorProjections cM nullVecs (LA.ident nNullVecs)
  where
    nProbs = S.size $ BRK.elements @k
    cM = DED.mMatrix nProbs $ DMS.msStencils ms
    nullVecs = nullSpaceVectors (DMS.msNumCategories ms) (DMS.msStencils ms)
    nNullVecs = fst $ LA.size nullVecs
-}

-- given list in a order, produce list in b order
permuteList :: forall a b c . DMS.IsomorphicKeys a b -> [c] -> [c]
permuteList ik cs = case ik of
  DMS.IsomorphicKeys abF _ -> fmap snd $ sortOn (abF . fst) $ zip (S.toList $ BRK.elements @a) cs

permutationMatrix :: forall a b . DMS.IsomorphicKeys a b -> LA.Matrix Double
permutationMatrix ik = case ik of
  DMS.IsomorphicKeys abF _ -> let aIndexL = zip (S.toList $ BRK.elements @a) [(0 :: Int)..]
                                  bIndexL = sortOn (abF . fst) aIndexL
                                  mAssoc = fmap (,1) $ zip [0..] $ fmap snd bIndexL
                                  n = length aIndexL
                              in LA.assoc (n, n) 0 mAssoc

mapNullVectorProjections :: DMS.IsomorphicKeys a b -> NullVectorProjections a -> NullVectorProjections b
mapNullVectorProjections ikab nva = case ikab of
  (DMS.IsomorphicKeys _ _) -> case nva of
    (NullVectorProjections nsp) -> let permM = permutationMatrix ikab
                                       newNSP = DNS.mkNullSpacePartition (DNS.knownM nsp LA.<> LA.tr permM) (DNS.unknownM nsp LA.<> LA.tr permM)
                                   in NullVectorProjections newNSP

    (PCAWNullVectorProjections nsp sV rM) -> let permM = permutationMatrix ikab
                                                 newNSP = DNS.mkNullSpacePartition (DNS.knownM nsp LA.<> LA.tr permM) (DNS.unknownM nsp LA.<> LA.tr permM)
                                               in PCAWNullVectorProjections newNSP sV rM

applyNSPWeights :: ProjectionsToDiff k -> LA.Vector LA.R -> LA.Vector LA.R -> LA.Vector LA.R
applyNSPWeights (RawDiff nvps) projWs pV = pV + projToFull nvps projWs
applyNSPWeights (AvgdDiff aM nvps) projWs pV = pV + aM LA.#> projToFull nvps projWs

type ObjectiveF = forall k . NullVectorProjections k -> LA.Vector Double -> LA.Vector Double -> LA.Vector Double -> (Double, LA.Vector Double)
type OptimalOnSimplexF r =  forall k . ProjectionsToDiff k -> VS.Vector Double -> VS.Vector Double -> K.Sem r (VS.Vector Double)

viaOptimalWeights :: K.KnitEffects r => ObjectiveF -> Double -> OptimalOnSimplexF r
viaOptimalWeights objF pEps ptd projWs prodV = do
  let n = VS.sum prodV
      pV = VS.map (/ n) prodV
  ows <- DED.mapPE $ optimalWeights objF pEps (nullVectorProjections ptd) projWs pV
  pure $ VS.map (* n) $ applyNSPWeights ptd ows pV

-- alpha - alpha'
euclideanNS :: ObjectiveF
euclideanNS _nvps projWs _ v =
  let x = (v - projWs)
  in (VS.sum $ VS.map (^ (2 :: Int)) x, 2 * x)


-- B(B^TB)^-1 (alpha - alpha')
euclideanFull :: ObjectiveF
euclideanFull nvps projWs _ v =
  let d = v - projWs
      a = projToFullM nvps
      x = a LA.#> d
  in (VS.sum $ VS.map (^ (2 :: Int)) x, LA.tr a LA.#> (2 * x))

euclideanWeighted :: (Double -> Double) -> ObjectiveF
euclideanWeighted g nvps projWs pV v =
  let d = v - projWs
      f y = if y < 1e-12 then 0 else g y
      invP = VS.map f pV
      a = LA.diag invP LA.<> projToFullM nvps
      x = a LA.#> d
  in (VS.sum $ VS.map (^ (2 :: Int)) x, LA.tr a LA.#> (2 * x))

klDiv :: ObjectiveF
klDiv nvps projWs pV v =
  let nV = pV + projToFull nvps projWs
      mV = pV + projToFull nvps v
  in (DED.klDiv nV mV, fullToProj nvps (DED.klGrad nV mV))

optimalWeights :: DED.EnrichDataEffects r
               => ObjectiveF
               -> Double
               -> NullVectorProjections k
               -> LA.Vector LA.R
               -> LA.Vector LA.R
               -> K.Sem r (LA.Vector LA.R)
optimalWeights objectiveF pEps nvps projWs pV = do
--  K.logLE K.Info $ "optimalWeights: pV = " <> DED.prettyVector pV
--  K.logLE K.Info $ "optimalWeights: Initial pV + nsWs <.> nVs = " <> DED.prettyVector (pV + projToFull nvps projWs)
  let n = VS.length projWs
--      prToFull = projToFull nvps
--      scaleGradM = fullToProjM nvps LA.<> LA.tr (fullToProjM nvps)
      objD = objectiveF nvps projWs pV
--      obj v = fst $ objD v
      constraintData =  L.zip (VS.toList pV) (LA.toRows $ projToFullM nvps)
      constraintLB :: (Double, LA.Vector LA.R)-> LA.Vector LA.R -> (Double, LA.Vector LA.R)
      constraintLB (p, projToNullC) v = (negate (p + v `LA.dot` projToNullC), negate projToNullC)
      constraintLBs = fmap constraintLB constraintData
      nlConstraintsD = fmap (\cf -> NLOPT.InequalityConstraint (NLOPT.Scalar cf) 1e-6) $ constraintLBs
--      nlConstraints = fmap (\cf -> NLOPT.InequalityConstraint (NLOPT.Scalar $ \v -> fst $ cf v) 1e-6) $ constraintLBs
{-      constraintV v = negate $ pV + projToFull nvps v
      constraintG = negate $ projToFullM nvps
      vecConstraintD v _ = (constraintV v , constraintG)
      nlvConstraintsD = [NLOPT.InequalityConstraint (NLOPT.Vector (fromIntegral $ VS.length pV) vecConstraintD) 1e-6]
-}
      -- we don't need upper bounds because those are enforced by a constraint.
      -- sum(pV + projToFullM w) = 1. And since all (pV + projToFull w)_i >= 0
      -- we know all (pV + projToFull w)_i <= 1
      maxIters = 1000
      absTol = pEps
      -- dV_j = sum_k B_{jk} d\alpha_k
      -- Assume d\alpha_k are iid: Var(dV_j) ~ Var(d\alpha) \sum_k B_{jk}^2
      -- d\alpha are uniform [-\eps, \eps] so Var(d\alpha) ~ \eps / \sqrt{3}
      -- \eps_\alpha = \eps/\sqrt 3 * \sum_k B_{jk}^2
      sdSumPtF = VS.fromList $ fmap (sqrt . VS.sum . VS.map (\x -> x * x)) . LA.toColumns $ projToFullM nvps
      absTolV = VS.map (\x -> absTol / (sqrt 3 * x)) sdSumPtF -- VS.fromList $ L.replicate n absTol
--  K.logLE K.Info $ "absTolV = " <> DED.prettyVector absTolV
  let nlStop = NLOPT.ParameterAbsoluteTolerance absTolV :| [NLOPT.MaximumEvaluations maxIters]
      nlAlgo = NLOPT.SLSQP objD [] nlConstraintsD [] --[nlSumToOneD]
--      nlAlgo = NLOPT.MMA objD nlConstraintsD
--      nlAlgo = NLOPT.COBYLA obj [] nlConstraints [] Nothing
      nlProblem =  NLOPT.LocalProblem (fromIntegral n) nlStop nlAlgo
      nlSol = NLOPT.minimizeLocal nlProblem projWs
  case nlSol of
    Left result -> PE.throw $ DED.TableMatchingException  $ "minConstrained: NLOPT solver failed: " <> show result
    Right solution -> case NLOPT.solutionResult solution of
      NLOPT.MAXEVAL_REACHED -> PE.throw $ DED.TableMatchingException $ "minConstrained: NLOPT Solver hit max evaluations (" <> show maxIters <> ")."
      NLOPT.MAXTIME_REACHED -> PE.throw $ DED.TableMatchingException $ "minConstrained: NLOPT Solver hit max time."
      _ -> do
        let oWs = NLOPT.solutionParams solution
--        K.logLE K.Info $ "solution=" <> DED.prettyVector oWs
--        K.logLE K.Info $ "Solution: pV + oWs <.> nVs = " <> DED.prettyVector (pV + projToFull nvps oWs)
        pure oWs


weightMapA :: (Double -> Double) -> LA.Vector Double -> LA.Matrix Double -> (LA.Matrix Double, [Int])
weightMapA g pV a =
  let f y = if y < 1e-12 then 0 else g y
      invP = VS.map f pV
      (_, nonZeroIs) = first (fmap fst) $ second (fmap fst) $ List.partition ((== 0) . snd) $ zip [0..] (VS.toList invP)
      aNZ = AS.subMatrixRL nonZeroIs (LA.diag invP LA.<> a)
  in (aNZ, nonZeroIs)

optimalWeightsAS :: DED.EnrichDataEffects r
                 => AS.ActiveSetConfiguration
                 -> Maybe (Double -> Double)
                 -> Maybe AS.LSI_E
                 -> NullVectorProjections k
                 -> LA.Vector LA.R
                 -> LA.Vector LA.R
                 -> K.Sem r (LA.Vector LA.R)
optimalWeightsAS asConfig mf mLSIE nvps projWs pV = do
  -- convert to correct form for AS solver
  let a = projToFullM nvps
  when (LA.cols a /= LA.size projWs)
    $ PE.throw $ DED.TableMatchingException
    $ "optimalWeightAS: cols(B^\\dagger) = " <> show (LA.cols a) <> " != " <> show (LA.size pV) <> " = length(pV)"
  let (aNZ, nonZeroIs) = maybe (a, [0..(LA.rows a - 1)]) (\f -> weightMapA f pV a) mf
      lsiE = fromMaybe (AS.Original aNZ) mLSIE
--      removeZeros =  VS.fromList . fmap snd . filter (not . (`elem` zeroIs) . fst) . zip [0..] . VS.toList
  let bNZ = aNZ LA.#> projWs
      pVNZ = AS.subVectorL nonZeroIs pV --removeZeros pV
      ic = AS.MatrixLower aNZ (negate pVNZ)
  (resE, _) <- AS.optimalLSI (K.logLE K.Info) asConfig lsiE bNZ ic
  case resE of
    Left err -> PE.throw $ DED.TableMatchingException $ "ActiveSet.findOptimal: " <> err
    Right ows -> pure ows

viaOptimalWeightsAS :: K.KnitEffects r => Maybe (Double -> Double) -> Maybe AS.LSI_E -> OptimalOnSimplexF r
viaOptimalWeightsAS mf mLSIE ptd projWs prodV = do
  let n = VS.sum prodV
      pV = VS.map (/ n) prodV
  ows <- DED.mapPE $ optimalWeightsAS AS.defaultActiveSetConfig mf mLSIE (nullVectorProjections ptd) projWs pV
  pure $ VS.map (* n) $ applyNSPWeights ptd ows pV


viaNearestOnSimplex :: OptimalOnSimplexF r
viaNearestOnSimplex nvps projWs prodV = do
  let n = VS.sum prodV
  pure $ VS.map (* n) $ projectToSimplex $ applyNSPWeights nvps projWs (VS.map (/ n) prodV)

viaNearestOnSimplexPlus :: DED.EnrichDataEffects r => OptimalOnSimplexF r
viaNearestOnSimplexPlus nvps projWs prodV = do
  let n = VS.sum prodV
      tgtV = projectToSimplex $ applyNSPWeights nvps projWs (VS.map (/ n) prodV)
  VS.map (* n) <$> optimalVector (nullVectorProjections nvps) prodV tgtV

-- after Chen & Ye: https://arxiv.org/pdf/1101.6081
projectToSimplex :: VS.Vector Double -> VS.Vector Double
projectToSimplex y = VS.fromList $ fmap (\x -> max 0 (x - tHat)) yL
  where
    yL = VS.toList y
    n = VS.length y
    sY = sort yL
    t i = (FL.fold FL.sum (L.drop i sY) - 1) / realToFrac (n - i)
    tHat = go (n - 1)
    go 0 = t 0
    go k = let tk = t k in if tk > sY L.!! k then tk else go (k - 1)

applyNSPWeightsO :: DED.EnrichDataEffects r => ObjectiveF -> Double -> ProjectionsToDiff k -> LA.Vector LA.R -> LA.Vector LA.R -> K.Sem r (LA.Vector LA.R)
applyNSPWeightsO objF pEps ptd nsWs pV = f <$> optimalWeights objF pEps (nullVectorProjections ptd) nsWs pV
  where f oWs = applyNSPWeights ptd oWs pV

-- this aggregates over cells with the same given key
labeledRowsToVecFld :: (Ord k, BRK.FiniteSet k, VS.Storable x, Num x, Monoid a) => Lens' a x -> (row -> k) -> (row -> a) -> FL.Fold row (VS.Vector x)
labeledRowsToVecFld wgtLens key dat
  = VS.fromList . fmap getSum . M.elems <$> (FL.premap (\r -> (key r, Sum $ view wgtLens $ dat r)) DMS.zeroFillSummedMapFld)
{-# SPECIALIZE labeledRowsToVecFld :: (Ord k, BRK.FiniteSet k, Monoid a) => Lens' a Double -> (row -> k) -> (row -> a) -> FL.Fold row (VS.Vector Double) #-}

labeledRowsToListFld :: (Ord k, BRK.FiniteSet k, Monoid w) => (row -> k) -> (row -> w) -> FL.Fold row [w]
labeledRowsToListFld key dat = M.elems <$> (FL.premap (\r -> (key r, dat r)) DMS.zeroFillSummedMapFld)


labeledRowsToNormalizedTableMapFld :: forall a outerK row w . (Ord a, BRK.FiniteSet a, Ord outerK, BRK.FiniteSet outerK, Monoid w)
                                   => Lens' w Double -> (row -> outerK) -> (row -> a) -> (row -> w) -> FL.Fold row (Map outerK (Map a w))
labeledRowsToNormalizedTableMapFld wgtLens outerKey innerKey nF = f <$> labeledRowsToKeyedListFld keyF nF where
  keyF row = (outerKey row, innerKey row)
  f :: [((outerK, a), w)] -> Map outerK (Map a w)
  f = FL.fold (DMS.normalizedTableMapFld wgtLens)

labeledRowsToTableMapFld :: forall a outerK row w . (Ord a, BRK.FiniteSet a, Ord outerK, BRK.FiniteSet outerK, Monoid w)
                         => (row -> outerK) -> (row -> a) -> (row -> w) -> FL.Fold row (Map outerK (Map a w))
labeledRowsToTableMapFld outerKey innerKey nF = f <$> labeledRowsToKeyedListFld keyF nF where
  keyF row = (outerKey row, innerKey row)
  f :: [((outerK, a), w)] -> Map outerK (Map a w)
  f = FL.fold DMS.tableMapFld

labeledRowsToKeyedListFld :: (Ord k, BRK.FiniteSet k, Monoid w) => (row -> k) -> (row -> w) -> FL.Fold row [(k, w)]
labeledRowsToKeyedListFld key dat =  M.toList <$> FL.premap (\r -> (key r, dat r)) DMS.zeroFillSummedMapFld
  -- fmap (zip (S.toList BRK.elements) . VS.toList) $ labeledRowsToVecFld key dat

labeledRowsToVec :: (Ord k, BRK.FiniteSet k, VS.Storable x, Num x, Foldable f, Monoid a)
                 => Lens' a x -> (row -> k) -> (row -> a) -> f row -> VS.Vector x
labeledRowsToVec wl key dat = FL.fold (labeledRowsToVecFld wl key dat)

labeledRowsToList :: (Ord k, BRK.FiniteSet k, Foldable f, Monoid w) => (row -> k) -> (row -> w) -> f row -> [w]
labeledRowsToList key dat = FL.fold (labeledRowsToListFld key dat)

productFromJointFld :: DMS.MarginalStructure w k -> (row -> k) -> (row -> w) -> FL.Fold row [(k, w)]
productFromJointFld ms keyF datF = case ms of
  DMS.MarginalStructure _ ptFld -> FL.fold ptFld . M.toList <$> FL.premap (\r -> (keyF r, datF r)) DMS.zeroFillSummedMapFld

productVecFromJointFld :: DMS.MarginalStructure w k -> Lens' w Double -> (row -> k) -> (row -> w) -> FL.Fold row (VS.Vector Double)
productVecFromJointFld ms wl keyF datF = fmap (VS.fromList . fmap (view wl . snd)) $ productFromJointFld ms keyF datF

diffProjectionsFromJointKeyedList :: forall k w . DMS.MarginalStructure w k
                                -> Lens' w Double
                                -> (VS.Vector Double -> VS.Vector Double)
                                -> [(k, w)]
                                -> VS.Vector Double
diffProjectionsFromJointKeyedList ms wl projDiff kl = projDiff (wgtVec kl - prodWeights kl)
  where
    wgtVec = VS.fromList . fmap (view wl . snd)
    prodWeights :: [(k, w)] -> VS.Vector Double
    prodWeights x = case ms of
      DMS.MarginalStructure _ ptFld -> wgtVec $ FL.fold ptFld x

diffProjectionsFromJointFld :: forall row k w . (BRK.FiniteSet k, Ord k, Monoid w)
                               => DMS.MarginalStructure w k
                               -> Lens' w Double
                               -> (VS.Vector Double -> VS.Vector Double)
                               -> (row -> k)
                               -> (row -> w)
                               -> FL.Fold row (VS.Vector Double)
diffProjectionsFromJointFld ms wl projDiff keyF datF = fmap (diffProjectionsFromJointKeyedList ms wl projDiff . M.toList)
                                                       $ FL.premap (\r -> (keyF r, datF r)) (DMS.normalizeAndFillMapFld wl)


sumLens :: Lens' (Sum x) x
sumLens = lens getSum (\_sx x -> Sum x)

applyNSPWeightsFld :: forall outerKs ks count rs r .
                      ( DED.EnrichDataEffects r
                      , V.KnownField count
                      , outerKs V.++ (ks V.++ '[count]) ~ (outerKs V.++ ks) V.++ '[count]
                      , ks F.⊆ (ks V.++ '[count])
                      , Integral (V.Snd count)
                      , F.ElemOf (ks V.++ '[count]) count
                      , Ord (F.Record ks)
                      , BRK.FiniteSet (F.Record ks)
                      , Ord (F.Record outerKs)
                      , outerKs F.⊆ rs
                      , (ks V.++ '[count]) F.⊆ rs
                      , FSI.RecVec (outerKs V.++ (ks V.++ '[count]))
                      )
                   => (F.Record outerKs -> Maybe Text)
                   -> ProjectionsToDiff (F.Record ks) -- ?
                   -> (F.Record outerKs -> FL.FoldM (K.Sem r) (F.Record (ks V.++ '[count])) (LA.Vector Double))
                   -> FL.FoldM (K.Sem r) (F.Record rs) (F.FrameRec (outerKs V.++ ks V.++ '[count]))
applyNSPWeightsFld logM ptd nsWsFldF =
  let keysL = S.toList $ BRK.elements @(F.Record ks) -- these are the same for each outerK
      precomputeFld :: F.Record outerKs -> FL.FoldM (K.Sem r) (F.Record (ks V.++ '[count])) (LA.Vector LA.R, Double, LA.Vector LA.R)
      precomputeFld ok =
        let sumFld = FL.generalize $ FL.premap (realToFrac . F.rgetField @count) FL.sum
            vecFld = FL.generalize $ labeledRowsToVecFld sumLens (F.rcast @ks) (Sum . realToFrac . F.rgetField @count)
        in (,,) <$> nsWsFldF ok <*> sumFld <*> vecFld
      compute :: F.Record outerKs -> (LA.Vector LA.R, Double, LA.Vector LA.R) -> K.Sem r [F.Record (outerKs V.++ ks V.++ '[count])]
      compute ok (nsWs, vSum, v) = do
        maybe (pure ()) (\msg -> K.logLE K.Info $ msg <> " nsWs=" <> DED.prettyVector nsWs) $ logM ok
        let optimalV = VS.map (* vSum) $ projectToSimplex $ applyNSPWeights ptd nsWs $ VS.map (/ vSum) v
--        optimalV <- fmap (VS.map (* vSum)) $ applyNSPWeightsO nvps nsWs $ VS.map (/ vSum) v
        pure $ zipWith (\k c -> ok F.<+> k F.<+> FT.recordSingleton @count (round c)) keysL (VS.toList optimalV)
      innerFld ok = FMR.postMapM (compute ok) (precomputeFld ok)
  in  FMR.concatFoldM
        $ FMR.mapReduceFoldM
        (FMR.generalizeUnpack $ FMR.noUnpack)
        (FMR.generalizeAssign $ FMR.assignKeysAndData @outerKs @(ks V.++ '[count]))
        (FMR.ReduceFoldM $ \k -> F.toFrame <$> innerFld k)


applyNSPWeightsFldG :: forall outerK k w r .
                      ( DED.EnrichDataEffects r
                      , Monoid w
                      , BRK.FiniteSet k
                      , Ord k
                      , Ord outerK
                      )
                    => Lens' w Double
                    -> (Double -> w -> w) -- not sure if I can make this update follow lens laws
                    -> (outerK -> Maybe Text)
                    -> ProjectionsToDiff k
                    -> (outerK -> FL.FoldM (K.Sem r) (k, w) (LA.Vector Double))
                    -> FL.FoldM (K.Sem r) (outerK, k, w) [(outerK, k, w)]
applyNSPWeightsFldG wl updateW logM ptd nsWsFldF =
  let keysL = S.toList $ BRK.elements @k -- these are the same for each outerK
      okF (ok, _, _) = ok
      kwF (_, k, w) = (k, w)
      precomputeFld :: outerK -> FL.FoldM (K.Sem r) (k, w) (LA.Vector Double, Double, [w])
      precomputeFld ok =
        let sumFld = FL.generalize $ FL.premap (view wl . snd) FL.sum
            vecFld = FL.generalize $ labeledRowsToListFld fst snd
        in (,,) <$> nsWsFldF ok <*> sumFld <*> vecFld
      compute :: outerK -> (LA.Vector LA.R, Double, [w]) -> K.Sem r [(outerK, k, w)]
      compute ok (nsWs, vSum, ws) = do
        maybe (pure ()) (\msg -> K.logLE K.Info $ msg <> " nsWs=" <> DED.prettyVector nsWs) $ logM ok
        let optimalV = VS.map (* vSum) $ projectToSimplex $ applyNSPWeights ptd nsWs $ VS.fromList $ fmap ((/ vSum) . view wl) ws
--        optimalV <- fmap (VS.map (* vSum)) $ applyNSPWeightsO nvps nsWs $ VS.map (/ vSum) v
        pure $ List.zipWith3 (\k c w -> (ok, k, updateW c w)) keysL (VS.toList optimalV) ws
      innerFldM ok = FMR.postMapM (compute ok) (precomputeFld ok)
  in  MR.concatFoldM
        $ MR.mapReduceFoldM
        (MR.generalizeUnpack MR.noUnpack)
        (MR.generalizeAssign $ MR.assign okF kwF)
        (FMR.ReduceFoldM innerFldM)


nullSpaceProjections :: (Ord k, Ord outerK, Show outerK, DED.EnrichDataEffects r, Foldable f)
                            => LA.Matrix LA.R
                            -> (row -> outerK)
                            -> (row -> k)
                            -> (row -> Double)
                            -> f row
                            -> f row
                            -> K.Sem r (Map outerK (Double, LA.Vector LA.R))
nullSpaceProjections nullVs outerKey key wgt actuals products = do
  let normalizedVec' m s = VS.fromList $ (/ s)  <$> M.elems m
      mapData m s = (s, normalizedVec' m s)
      toDatFld = mapData <$> FL.map <*> FL.premap snd FL.sum
      toMapFld = FL.premap (\r -> (outerKey r, (key r, wgt r))) $ FL.foldByKeyMap toDatFld
      actualM = FL.fold toMapFld actuals
      prodM = FL.fold toMapFld products
      whenMatchedF _ (na, aV) (_, pV) = pure (na, nullVs LA.#> (aV - pV))
      whenMissingAF outerK _ = PE.throw $ DED.TableMatchingException $ "averageNullSpaceProjections: Missing actuals for outerKey=" <> show outerK
      whenMissingPF outerK _ = PE.throw $ DED.TableMatchingException $ "averageNullSpaceProjections: Missing product for outerKey=" <> show outerK
  MM.mergeA (MM.traverseMissing whenMissingAF) (MM.traverseMissing whenMissingPF) (MM.zipWithAMatched whenMatchedF) actualM prodM

avgNullSpaceProjections :: (Double -> Double) -> Map a (Double, LA.Vector LA.R) -> LA.Vector LA.R
avgNullSpaceProjections popWeightF m = VS.map (/ totalWeight) . VS.fromList . fmap (FL.fold FL.mean . VS.toList) . LA.toColumns . LA.fromRows . fmap weight $ M.elems m
  where
    totalWeight = FL.fold (FL.premap (popWeightF . fst) FL.sum) m / realToFrac (M.size m)
    weight (n, v) = VS.map (* popWeightF n) v

diffCovarianceFldMS :: forall outerK k row w .
                       (Ord outerK, BRK.FiniteSet k)
                    => Lens' w Double
                    -> (row -> outerK)
                    -> (row -> k)
                    -> (row -> w)
                    -> DMS.MarginalStructure w k
                    -> Maybe (DNS.CatsAndKnowns k)
                    -> Maybe (FL.Fold row (LA.Vector Double, LA.Herm Double))
diffCovarianceFldMS wl outerKey catKey dat ms cmM = do
  nVs <- nvpProj <$> nullVecsMS ms cmM
  pure $ case ms of
           (DMS.MarginalStructure subsets ptFld) -> diffCovarianceFld wl outerKey catKey dat
                                                    nVs
                                                    (fmap snd . FL.fold ptFld)

diffCovarianceFld :: forall outerK k row w .
                     (Ord outerK, Ord k, BRK.FiniteSet k, Monoid w)
                  => Lens' w Double
                  -> (row -> outerK)
                  -> (row -> k)
                  -> (row -> w)
                  -> LA.Matrix LA.R
                  -> ([(k, w)] -> [w])
                  -> FL.Fold row (LA.Vector LA.R, LA.Herm LA.R)
diffCovarianceFld wl outerKey catKey dat nullVecs' prodTableF = LA.meanCov . LA.fromRows <$> vecsFld
  where
    allKs :: Set k = BRK.elements
--    n = S.size allKs
--    (_, nullVecs') = nullSpaceDecomposition n sts
    wgtVec = VS.fromList . fmap (view wl)
    pcF :: [w] -> VS.Vector Double
    pcF = wgtVec . prodTableF . zip (S.toList allKs)
    projections ws = let ws' = DMS.normalize wl ws in nullVecs' LA.#> (wgtVec ws' - pcF ws')
    innerFld = projections <$> labeledRowsToListFld fst snd
    vecsFld = MR.mapReduceFold
              MR.noUnpack
              (MR.assign outerKey (\row -> (catKey row, dat row)))
              (MR.foldAndLabel innerFld (\_ v -> v))

nullSpaceDecompositionMS :: forall k w . BRK.FiniteSet k => DMS.MarginalStructure w k -> DNS.NullSpacePartition
nullSpaceDecompositionMS = \cases
  (DMS.MarginalStructure subsets _) -> nullSpaceDecompositionSubsets @k subsets

nullSpaceDecompositionSubsets :: forall k . (BRK.FiniteSet k, Ord k) => [Set k] -> DNS.NullSpacePartition
nullSpaceDecompositionSubsets subsets = DNS.nullSpacePartitionSVD (S.size $ BRK.elements @k) $ fmap DMS.subsetToStencil subsets

uncorrelatedNullVecsMS :: forall k w . DMS.MarginalStructure w k -> Maybe (DNS.CatsAndKnowns k) -> LA.Herm LA.R -> Maybe (NullVectorProjections k)
uncorrelatedNullVecsMS ms cmM cov = do
  plainNVPs <- nullVecsMS ms cmM
  pure $ case ms of
    (DMS.MarginalStructure _ _) -> uncorrelatedNullVecs (nvpConstraints plainNVPs) (nvpProj plainNVPs) cov

uncorrelatedNullVecs :: (Ord k, BRK.FiniteSet k)
                     => LA.Matrix LA.R
                     -> LA.Matrix LA.R
                     -> LA.Herm LA.R
                     -> NullVectorProjections k
uncorrelatedNullVecs cMatrix nVs cov = PCAWNullVectorProjections (DNS.mkNullSpacePartition cMatrix nVs) eigVals eigVecs
  where
    (eigVals, eigVecs) = LA.eigSH cov

nullVecsMS :: forall k w . DMS.MarginalStructure w k  -> Maybe (DNS.CatsAndKnowns k) -> Maybe (NullVectorProjections k)
nullVecsMS ms cmM = case ms of
   (DMS.MarginalStructure _ _) -> NullVectorProjections  <$> cnM
     where
       cnM = case cmM of
         Nothing ->  Just $ nullSpaceDecompositionMS ms
         Just catsAndMarginals ->  DNS.nullSpacePartitionCM catsAndMarginals

----
optimalVector :: DED.EnrichDataEffects r
               => NullVectorProjections k
               -> LA.Vector LA.R
               -> LA.Vector LA.R
               -> K.Sem r (LA.Vector LA.R)
optimalVector nvps pV tgtV = do
--  K.logLE K.Info $ "Initial: V=" <> DED.prettyVector pV
--  K.logLE K.Info $ "Tgt: V=" <> DED.prettyVector tgtV
  let n = VS.length pV
--      prToFull = projToFull nvps
  --      scaleGradM = fullToProjM nvps LA.<> LA.tr (fullToProjM nvps)
--      objD v = (DED.klDivP v tgtV, negate $ DED.klGradP' v tgtV)
      objD v = let d =  VS.zipWith (-) v tgtV in (LA.norm_2 d, 2 * d)
      cM = nvpConstraints nvps
      nVs = fullToProjM nvps
--  K.logLE K.Info $ "Constraints: " <> show (LA.size cM)
--  K.logLE K.Info $ "Null-Space basis" <> show (LA.size nVs)
--  let (u, _, _) = LA.compactSVD $ LA.tr cM
  let cRHS = cM LA.#> pV
      constraintData = L.zip (VS.toList cRHS) (LA.toRows cM)
      constraintF (cVal, cRow) v = ((cRow `LA.dot` v) - cVal, cRow)
      constraintFs = fmap constraintF constraintData
      nlConstraintsD = fmap (\cf -> NLOPT.EqualityConstraint (NLOPT.Scalar cf) 1e-6) constraintFs
      lowerBounds = NLOPT.LowerBounds $ VS.replicate n 0
      upperBounds = NLOPT.UpperBounds $ VS.replicate n 1
--      nlConstraints = fmap (\cf -> NLOPT.InequalityConstraint (NLOPT.Scalar $ fst . cf) 1e-5) constraintFs
      maxIters = 1000
      absTol = 1e-6
      absTolV = VS.fromList $ L.replicate n absTol
      nlStop = NLOPT.ParameterAbsoluteTolerance absTolV :| [NLOPT.MaximumEvaluations maxIters]
      nlAlgo = NLOPT.SLSQP objD [lowerBounds, upperBounds] [] nlConstraintsD
--      nlAlgo = NLOPT.MMA objD nlConstraintsD
      nlProblem =  NLOPT.LocalProblem (fromIntegral n) nlStop nlAlgo
      nlSol = NLOPT.minimizeLocal nlProblem pV
  case nlSol of
    Left result -> PE.throw $ DED.TableMatchingException  $ "minConstrained: NLOPT solver failed: " <> show result
    Right solution -> case NLOPT.solutionResult solution of
      NLOPT.MAXEVAL_REACHED -> PE.throw $ DED.TableMatchingException $ "minConstrained: NLOPT Solver hit max evaluations (" <> show maxIters <> ")."
      NLOPT.MAXTIME_REACHED -> PE.throw $ DED.TableMatchingException $ "minConstrained: NLOPT Solver hit max time."
      _ -> do
        let oWs = NLOPT.solutionParams solution
--        K.logLE K.Info $ "solution=" <> DED.prettyVector oWs
--        K.logLE K.Info $ "Solution: pV + oWs <.> nVs = " <> DED.prettyVector (pV + projToFull nvps oWs)
        pure oWs
