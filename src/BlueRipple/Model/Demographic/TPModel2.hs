{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE UnicodeSyntax #-}
{-# LANGUAGE StandaloneDeriving #-}

module BlueRipple.Model.Demographic.TPModel2
  (
    module BlueRipple.Model.Demographic.TPModel2
  )
where

import Relude.Extra (traverseToSnd)

import qualified BlueRipple.Data.CachingCore as BRCC
--import qualified BlueRipple.Data.CachingCore as BRK

import qualified BlueRipple.Model.Demographic.DataPrep as DDP
import qualified BlueRipple.Model.Demographic.MarginalStructure as DMS
import qualified BlueRipple.Model.Demographic.TableProducts as DTP
import qualified BlueRipple.Model.StanTools as MST

import qualified BlueRipple.Data.Keyed as BRK
import qualified BlueRipple.Data.Types.Demographic as DT
import qualified BlueRipple.Data.Types.Geographic as GT
import qualified BlueRipple.Data.ACS_PUMS as ACS

import qualified Knit.Report as K

import qualified Control.MapReduce.Simple as MR

import qualified Control.Foldl as FL
import qualified Data.Map.Strict as M
import qualified Data.Set as S

import qualified Data.List as List
import qualified Frames as F
import qualified Frames.Serialize as FS
import qualified Numeric.LinearAlgebra as LA
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vinyl as V

import Control.Lens (view, Lens')
import GHC.TypeLits (Symbol)

import qualified Stan as S
import qualified Stan.BuildingBlocks as SBB (rowLength)
import Stan (TypedList(..))
import Stan.Operators
--import qualified CmdStan as CS

import qualified Flat

productDistributionFld :: forall outerK k row w .
                          (Ord outerK)
                       => DMS.MarginalStructure w k
                       -> Lens' w Double
                       -> (row -> outerK)
                       -> (row -> k)
                       -> (row -> w)
                       -> FL.Fold row (Map outerK (VS.Vector Double))
productDistributionFld marginalStructure wl outerKey catKey datF = M.fromList <$> case marginalStructure of
  DMS.MarginalStructure _ ptFld -> MR.mapReduceFold
                                   MR.noUnpack
                                   (MR.assign outerKey id)
                                   (MR.foldAndLabel innerFld (,))
    where
      pcF =  VS.fromList . fmap (view wl . snd) . FL.fold ptFld
      innerFld = DTP.normalizedVec . pcF <$> DTP.labeledRowsToKeyedListFld catKey datF



-- produce the projections of the difference bewteen the distirbution of
-- probability 1 at k and the product distribution at outerK
rowDiffProjections ::  forall outerK k row .
                       (Ord outerK, Show outerK, Ord k, BRK.FiniteSet k)
                   => DTP.NullVectorProjections k
                   -> Map outerK (VS.Vector Double) -- Product Distribution
                   -> (row -> outerK)
                   -> (row -> k)
                   -> row
                   -> Either Text (VS.Vector Double)
rowDiffProjections nvps pdMap outerKey catKey r = do
  let ok = outerKey r
      k = catKey r
      catDist = VS.fromList $ fmap getSum $ M.elems (M.singleton k (Sum 1) <> DMS.zeroMap)
  pd <- maybeToRight ("rowDiffProjections: " <> show ok <> " is missing from product structure map!") $ M.lookup ok pdMap
  pure $ DTP.fullToProj nvps (catDist - pd)

rowsWithProjectedDiffs :: (Traversable g
                          , Ord outerK
                          , Show outerK
                          , Ord k
                          , BRK.FiniteSet k)
                       =>  DTP.NullVectorProjections k
                       -> Map outerK (VS.Vector Double) -- Product Distribution
                       -> (F.Record rs -> outerK)
                       -> (F.Record rs -> k)
                       -> g (F.Record rs)
                       -> Either Text (g (ProjDataRow rs))
rowsWithProjectedDiffs nvps pdMap outerKey catKey =
  fmap (fmap $ fmap (\(r, v) -> ProjDataRow r v))
  $ traverse (traverseToSnd $ rowDiffProjections nvps pdMap outerKey catKey)

data ProjDataRow rs = ProjDataRow (F.Record rs) (VS.Vector Double)

projRowRec :: ProjDataRow rs -> F.Record rs
projRowRec (ProjDataRow r _) = r

projRowVec :: ProjDataRow rs -> VS.Vector Double
projRowVec (ProjDataRow _ v) = v

instance (V.RMap rs, FS.RecFlat rs) => Flat.Flat (ProjDataRow rs) where
  size (ProjDataRow r v) = Flat.size (FS.toS r, VS.toList v)
  encode (ProjDataRow r v) = Flat.encode (FS.toS r, VS.toList v)
  decode = fmap (\(sr, l) -> ProjDataRow (FS.fromS sr) (VS.fromList l)) Flat.decode

data ProjData rs = ProjData {pdNNullVecs :: Int, pdNPredictors :: Int, pdRows :: [ProjDataRow rs]} deriving stock (Generic)
deriving anyclass instance (Flat.Flat (ProjDataRow rs)) => Flat.Flat (ProjData rs)

modelIDT :: S.InputDataType S.ModelDataT (ProjData rs)
modelIDT = S.ModelData

type ProjDataRTT rs = S.RowTypeTag (ProjDataRow rs)

data SlopeIntercept = SlopeIntercept { siSlope :: Double, siIntercept :: Double} deriving stock (Show, Generic)

applySlopeIntercept :: SlopeIntercept -> Double -> Double
applySlopeIntercept (SlopeIntercept s i) x = i + s * x
{-# INLINEABLE applySlopeIntercept #-}

newtype ModelResult g k (pd :: Type -> Type) = ModelResult { unModelResult :: Map g (Map k [Double], pd [SlopeIntercept]) }
  deriving stock (Generic)

deriving stock instance (Show g, Show k, Show (b [SlopeIntercept])) => Show (ModelResult g k b)
deriving anyclass instance (Ord g, Flat.Flat g, Ord k, Flat.Flat k, Flat.Flat (b [SlopeIntercept])) => Flat.Flat (ModelResult g k b)

modelResultNVPs :: (Traversable pd, Applicative pd, Show g, Ord g, Show k, Ord k)
                => ModelResult g k pd
                -> (r -> g)
                -> (r -> k)
                -> (r -> pd Double)
                -> r -> Either Text (VS.Vector Double)
modelResultNVPs modelResult geoKey catKey pdF r = do
  let gk = geoKey r
      ck = catKey r
      pd = pdF r
  (gaM, pdSIs) <- maybeToRight ("modelResultNVPs: " <> show gk <> " not found in model result geo-alpha map!")
        $ M.lookup gk $ unModelResult modelResult
  alphaV <- maybeToRight ("modelResultNVPs: " <> show ck <> " not found in model result alpha map for " <> show gk <> "!")
            $ M.lookup ck gaM
  let pdSIL = sequenceA pdSIs
      applyTo si = applySlopeIntercept <$> si <*> pd
      betaV = VS.fromList $ fmap (getSum . foldMap Sum . applyTo) pdSIL
  pure $ VS.fromList alphaV + betaV

stateG :: S.GroupTypeTag Text
stateG = S.GroupTypeTag "State"

stateGroupBuilder :: forall f rs . (Foldable f, Typeable rs)
                  => (F.Record rs -> Text) -> f Text -> S.StanDataBuilderEff S.ModelDataT (ProjData rs) (ProjDataRTT rs)
stateGroupBuilder saF states = do
  projData <- S.addData "ProjectionData" (modelIDT @rs) (S.ToFoldable pdRows)
  S.addGroupIndexForData (modelIDT @rs) stateG projData $ S.makeIndexFromFoldable show (saF . projRowRec) states
  S.addGroupIntMapForData (modelIDT @rs) stateG projData $ S.dataToIntMapFromFoldable (saF . projRowRec) states
  pure projData

data ProjModelData r =
  ProjModelData
  {
    projDataTag :: ProjDataRTT r
  , nNullVecsE :: S.IntE
  , nAlphasE :: S.IntE
  , alphasE :: S.MatrixE
  , nPredictorsE :: S.IntE
  , predictorsE :: S.MatrixE
  , projectionsE :: S.MatrixE
  , countsE :: S.IntArrayE
  }

data AlphaModel = AlphaSimple | AlphaHierCentered | AlphaHierNonCentered deriving stock (Show)

alphaModelText :: AlphaModel -> Text
alphaModelText AlphaSimple = "AS"
alphaModelText AlphaHierCentered = "AHC"
alphaModelText AlphaHierNonCentered = "AHNC"

data Distribution = NormalDist -- | CauchyDist | StudentTDist

distributionText :: Distribution -> Text
distributionText NormalDist = "normal"
--distributionText CauchyDist = "cauchy"
--distributionText StudentTDist = "studentT"

data ModelConfig fullK alphaK pd where
  ModelConfig :: Traversable pd
              => { projVecs :: DTP.NullVectorProjections fullK
                 , standardizeNVs :: Bool
                 , alphaDMR :: S.DesignMatrixRow alphaK
                 , predDMR :: S.DesignMatrixRow (pd Double)
                 , alphaModel :: AlphaModel
                 , distribution :: Distribution
                 } -> ModelConfig fullK alphaK pd

modelNumNullVecs :: ModelConfig fullK alphaK md -> Int
modelNumNullVecs mc = fst $ LA.size $ DTP.nvpProj mc.projVecs

modelText :: ModelConfig fullK alphaK md -> Text
modelText mc = distributionText mc.distribution <> "_" <> mc.alphaDMR.dmName <> "_" <> mc.predDMR.dmName <> "_" <> alphaModelText mc.alphaModel

dataText :: ModelConfig fullK alphaK md -> Text
dataText mc = mc.alphaDMR.dmName <> "_" <> mc.predDMR.dmName <> "_NV" <> show (modelNumNullVecs mc)

projModelData :: forall pd alphaK fullK rs . (Typeable rs)
              =>  ModelConfig fullK alphaK pd
              -> (F.Record rs -> alphaK)
              -> (F.Record rs -> Int)
              -> (F.Record rs -> pd Double)
              -> ProjDataRTT rs
              -> S.StanModelBuilderEff (ProjData rs) () (ProjModelData rs)
projModelData mc catKey countF predF projData = do
--  projData <- S.dataSetTag @(ProjDataRow rs) S.ModelData "ProjectionData"
  let projMER :: S.MatrixRowFromData (ProjDataRow r) --(outerK, md Double, VS.Vector Double)
      projMER = S.MatrixRowFromData "nvp" Nothing (modelNumNullVecs mc) (\(ProjDataRow _ v) -> VU.convert v)
      -- convert is here because we want unboxed vectors for JSON but hmatix uses storable vectors for FFI
  (pmE, nNullVecsE') <- S.add2dMatrixData (modelIDT @rs) projData projMER Nothing Nothing
--  let nNullVecsE' = S.mrfdColumnsE projMER
  let (_, nAlphasE') = S.designMatrixColDimBinding mc.alphaDMR Nothing
  alphaDME <- if SBB.rowLength mc.alphaDMR > 0
              then S.addDesignMatrix (modelIDT @rs) projData (contramap (catKey . projRowRec) mc.alphaDMR) Nothing
              else pure $ S.namedE "ERROR" S.SMat -- this shouldn't show up in stan code at all
  let (_, nPredictorsE') = S.designMatrixColDimBinding mc.predDMR Nothing
  dmE <- if SBB.rowLength mc.predDMR > 0
         then S.addDesignMatrix (modelIDT @rs) projData (contramap (predF . projRowRec) mc.predDMR) Nothing
         else pure $ S.namedE "ERROR" S.SMat -- this shouldn't show up in stan code at all
  countsE' <- S.addCountData (modelIDT @rs) projData "count" (countF . projRowRec)
  pure $ ProjModelData projData nNullVecsE' nAlphasE' alphaDME nPredictorsE' dmE pmE countsE'

-- S states
-- K projections
-- C categories
-- D predictors
-- either an K row-vector or S x K matrix
--data Alpha0 = SimpleAlpha0 S.RVectorE | HierarchicalAlpha0 S.MatrixE
-- C x K matrix or array[S] of C x K matrix
data Alpha = SimpleAlpha (S.Parameter S.EMat) | HierarchicalAlpha (S.Parameter (S.EArray1 S.EMat))
-- D x K matrix or Nothing
newtype Theta = Theta (Maybe (S.Parameter S.EMat))
-- sigma is K row-vector
newtype Sigma = Sigma { unSigma :: S.Parameter S.ERVec }

data ProjModelParameters where
  NormalParameters :: Alpha -> Theta -> Sigma -> ProjModelParameters

paramTheta :: ProjModelParameters -> Theta
paramTheta (NormalParameters _ t _) = t

projModelAlpha :: ModelConfig fullK alphaK pd -> ProjModelData rs -> S.StanModelBuilderEff (ProjData rs) () Alpha
projModelAlpha mc pmd = do
  let nStatesE = S.groupSizeE stateG
      hierAlphaSpec = S.array1Spec nStatesE (S.matrixSpec pmd.nAlphasE pmd.nNullVecsE)
      hierAlphaNDS = S.NamedDeclSpec "alpha" hierAlphaSpec
      indexAK a k = S.slice0 k . S.slice0 a
      indexSAK s a k = S.slice0  s . indexAK a k
      loopSAK stmtsF =
        S.nestedLoops (S.vftSized "s" nStatesE :> S.vftSized "a" pmd.nAlphasE :> S.vftSized "k" pmd.nNullVecsE :> TNil)
        $ \(s :> a :> k :> TNil) -> S.grouped $ stmtsF s a k

--      diagPostMult m cv = S.functionE S.diagPostMultiply (m :> cv :> TNil)
--      rowsOf nRowsE rv = diagPostMult (S.functionE S.rep_matrix (S.realE 1 :> nRowsE :> S.functionE S.size (rv :> TNil) :> TNil)) (S.transposeE rv)
      hierAlphaPs = do
        muAlphaP <- S.iidMatrixP
                    (S.NamedDeclSpec "muAlpha" $ S.matrixSpec pmd.nAlphasE pmd.nNullVecsE)
                    [] TNil S.std_normal
        sigmaAlphaP <- S.iidMatrixP
                       (S.NamedDeclSpec "sigmaAlpha" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) $ S.matrixSpec pmd.nAlphasE pmd.nNullVecsE)
                       [] TNil S.std_normal
        pure (muAlphaP :> sigmaAlphaP :> TNil)

  case mc.alphaModel of
    AlphaSimple -> do
      fmap SimpleAlpha
        $  S.iidMatrixP
        (S.NamedDeclSpec "alpha" $ S.matrixSpec pmd.nAlphasE pmd.nNullVecsE)
        [] TNil S.std_normal
    AlphaHierCentered -> do
      alphaPs <- hierAlphaPs
      fmap HierarchicalAlpha
        $ S.addBuildParameter
        $ S.UntransformedP hierAlphaNDS [] alphaPs
        $ \(muAlphaE :> sigmaAlphaE :> TNil) m
          -> S.addStmt
             $ loopSAK $ \s a k -> [S.sample (indexSAK s a k m)  S.normalS (indexAK a k muAlphaE :> indexAK a k sigmaAlphaE :> TNil)]
    AlphaHierNonCentered -> do
      alphaPs <- hierAlphaPs
      let rawNDS = S.NamedDeclSpec (S.rawName $ S.declName hierAlphaNDS) hierAlphaSpec

      rawP <- S.addBuildParameter
              $ S.UntransformedP rawNDS [] TNil
              $ \_ m -> S.addStmt $ loopSAK $ \s a k -> [S.sample (indexSAK s a k m) S.std_normal TNil]
      fmap HierarchicalAlpha
        $ S.addBuildParameter
        $ S.simpleTransformedP hierAlphaNDS [] (rawP :> alphaPs) S.TransformedParametersBlock
        $ \(rmE :> muE :> sigmaE :> TNil) ->
            let inner pE s a k = [indexSAK s a k pE S.|=| (indexAK a k muE |+| (indexAK a k sigmaE |*| indexSAK s a k rmE))]
            in S.DeclCodeF $ S.addStmt . loopSAK . inner

projModelParameters :: ModelConfig fullK alphaK pd -> ProjModelData rs -> S.StanModelBuilderEff (ProjData rs) () ProjModelParameters
projModelParameters mc pmd = do
  let stdNormalDWA :: (S.TypeOneOf t [S.EReal, S.ECVec, S.ERVec], S.GenSType t) => S.DensityWithArgs t
      stdNormalDWA = S.DensityWithArgs S.std_normal TNil --(S.realE 0 :> S.realE 1 :> TNil)
      numPredictors = SBB.rowLength mc.predDMR
  theta <- if numPredictors > 0 then
             fmap (Theta . Just)
               $ S.iidMatrixP
               (S.NamedDeclSpec "theta" $ S.matrixSpec pmd.nPredictorsE pmd.nNullVecsE)
               [] TNil
               S.std_normal
           else pure $ Theta Nothing
  sigma <-  fmap Sigma
             $ S.simpleParameterWA
             (S.NamedDeclSpec "sigma" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) $ S.rowVectorSpec pmd.nNullVecsE)
             stdNormalDWA
  alpha <- projModelAlpha mc pmd
  pure $ NormalParameters alpha theta sigma

data RunConfig = RunConfig { rcIncludePPCheck :: Bool, rcIncludeLL :: Bool }

projModel :: Typeable rs
          => RunConfig
          -> (F.Record rs -> alphaK)
          -> (F.Record rs -> Int)
          -> (F.Record rs -> pd Double)
          -> ModelConfig fullK alphaK pd
          -> ProjDataRTT rs
          -> S.StanModelBuilderEff  (ProjData rs) () ()
projModel rc alphaKey countF predF mc projData = do
  mData <- projModelData mc alphaKey countF predF projData
  mParams <- projModelParameters mc mData
  let betaNDS = S.NamedDeclSpec "beta" $ S.matrixSpec mData.nPredictorsE mData.nNullVecsE
      nRowsE = S.dataSetSizeE mData.projDataTag
      loopNVs = S.loopSized mData.nNullVecsE "k" --S.for "k" (S.SpecificNumbered (S.intE 1) mData.nNullVecsE)
      pExpr = S.parameterExpr
  -- transformedData
  (predM, _centerF, _mBeta) <- case paramTheta mParams of
    Theta (Just thetaP) -> do
      (centeredPredictorsE, centerF) <- S.centerDataMatrix S.DMCenterOnly mData.predictorsE Nothing "DM"
      (dmQ, _, _, mBeta) <- S.thinQR centeredPredictorsE "DM" $ Just (pExpr thetaP, betaNDS)
      pure (dmQ, centerF, mBeta)
    Theta Nothing -> pure (S.namedE "ERROR" S.SMat, \_ x _ -> pure x, Nothing)
  (nvps, inverseF) <- case mc.standardizeNVs of
    True -> S.inBlock S.SBTransformedData $ S.addFromCodeWriter $ do
      let nvVecDS t = S.NamedDeclSpec t $ S.rowVectorSpec mData.nNullVecsE
      sds <- S.declareNW (nvVecDS "nvpSDs")
      stdNVPs <- S.declareNW (S.NamedDeclSpec "stdNVPs" $ S.matrixSpec nRowsE mData.nNullVecsE)
      S.addStmt
        $ loopNVs
        $ \k -> S.grouped [ (sds !! k) S.|=| S.sd (mData.projectionsE `S.atCol` k)
                          , stdNVPs `S.atCol` k S.|=| ((mData.projectionsE `S.atCol` k) |/| (sds !! k))]
      let inverse :: (t ~ S.BinaryResultT S.BMultiply S.EReal t) => S.IntE -> S.UExpr t -> S.UExpr t --S.UExpr (S.BinaryResultT S.BMultiply S.EReal t)
          inverse k psCol = sds !! k |*| psCol
      pure (stdNVPs, inverse)
    False -> pure (mData.projectionsE, const id)
  countsVec <- S.inBlock S.SBTransformedData
               $ S.addFromCodeWriter
               $ S.declareRHSNW (S.NamedDeclSpec "countsV" $ S.vectorSpec nRowsE)
               $ S.to_vector mData.countsE
  -- model
  let reIndexByState = S.indexE S.s0 (S.dataByGroupIndexE mData.projDataTag stateG)
      -- given alpha and theta return an nData x nNullVecs matrix
      muE :: Alpha -> Theta -> S.StanModelBuilderEff (ProjData r) () S.MatrixE
      muE a t = S.addFromCodeWriter $ do
        let mThetaE = case t of
              Theta x -> fmap (\y -> predM |*| (pExpr y)) $ x
            muSpec = S.NamedDeclSpec "mu" $ S.matrixSpec nRowsE mData.nNullVecsE
        case a of
          SimpleAlpha alphaP -> case mThetaE of
            Nothing -> S.declareRHSNW muSpec (alphasE mData |*| pExpr alphaP)
            Just x -> S.declareRHSNW muSpec (alphasE mData |*| pExpr alphaP |+| x)
          HierarchicalAlpha alphaP -> do
            mu <- S.declareNW muSpec
            S.addStmt $ S.loopSized nRowsE "n" $ \n ->
              (mu `S.atRow` n)
                `S.assign`
                (case mThetaE of
                    Nothing -> (mData.alphasE `S.atRow` n) |*| ((reIndexByState $ pExpr alphaP) !! n)
                    Just mt -> (mData.alphasE `S.atRow` n) |*| ((reIndexByState $ pExpr alphaP) !! n) |+| (mt `S.atRow` n)
                )
            pure mu

      sigmaE :: Sigma -> S.IntE -> S.VectorE
      sigmaE s k = S.rep_vector (pExpr (unSigma s) !! k) nRowsE

      ppF :: S.MatrixE
          -> Int
          -> ((S.IntE -> S.ExprList xs) -> S.IntE -> S.UExpr S.EReal)
          -> (S.MatrixE -> S.IntE -> S.CodeWriter (S.IntE -> S.ExprList xs))
          -> S.StanModelBuilderEff (ProjData r) () (S.ArrayE S.EReal)
      ppF muMat k rngF rngPSCW = S.generatePosteriorPrediction'
                                 mData.projDataTag
                                 (S.NamedDeclSpec ("predProj_" <> show k) $ S.array1Spec nRowsE S.realSpec)
                                 rngF
                                 (rngPSCW muMat (S.intE k))
                                 (\_ p -> inverseF (S.intE k) p)
--      eltTimes = S.binaryOpE (S.SElementWise S.SMultiply)
      (muMatBuilder, sampleStmtF, ppStmtF) = case mParams of
        NormalParameters a t s ->
          let ssF e muMat k = S.familySample S.normalDist e (muMat `S.atCol` k :> sigmaE s k :> TNil)
                --S.familySample S.countScaledNormalDist e (countsVec :> muMat `S.atCol` k :> sigmaE s k :> TNil)
              rF f nE = S.familyRNG S.countScaledNormalDist (f nE) --S.functionE S.normal_rng (f nE)
              rpF muMat k = pure $ \nE -> countsVec !! nE :> muMat `S.atCol` k !! nE :> sigmaE s k !! nE :> TNil
          in (muE a t, ssF, \muMat n -> ppF muMat n rF rpF)

  S.inBlock S.SBModel $ do
    muMat <- muMatBuilder
    S.addFromCodeWriter $ do
      let loopBody k = S.cwStmt_ $ S.addStmt $ sampleStmtF (nvps `S.atCol` k) muMat k
      S.addStmt $ loopNVs loopBody
  -- generated quantities
  when rc.rcIncludePPCheck $ do
    muMat <- S.inBlock S.SBGeneratedQuantities muMatBuilder
    forM_ [1..modelNumNullVecs mc] (ppStmtF muMat)
  pure ()


runProjModel :: forall (ksO :: [(Symbol, Type)]) ksM pd r .
                (K.KnitEffects r
                , BRCC.CacheEffects r
                , ksM F.⊆ DDP.ACSa5ByPUMAR
                , ksO F.⊆ DDP.ACSa5ByPUMAR
                , Typeable pd
                , Ord (F.Record ksO)
                , BRK.FiniteSet (F.Record ksO)
--                , Flat.Flat (pd [SlopeIntercept])
                )
             => Bool
             -> Maybe Int
             -> RunConfig
             -> ModelConfig (F.Record ksO) (F.Record ksM) pd
             -> DMS.MarginalStructure (Sum Double) (F.Record ksO)
             -> (F.Record DDP.ACSa5ByPUMAR -> pd Double)
             -> K.Sem r (K.ActionWithCacheTime r ())
runProjModel clearCaches thinM rc mc ms predF = do
  let cacheRoot = "model/demographic/nullVecProjModel/"
      cacheDirE = (if clearCaches then Left else Right) cacheRoot
      dataName = "projectionData_" <> dataText mc
  stanDir <- K.liftKnit MST.stanDir >>= K.knitMaybe "runModel: empty stanDir!" . BRCC.insureFinalSlash
  let runnerInputNames = S.RunnerInputNames
                         (stanDir <> "demographic/nullVecProj2")
                         (modelText mc)
                         (Just $ S.GQNames "pp" dataName) -- posterior prediction vars to wrap
                         dataName
      (srcWindow, cachedSrc) = ACS.acs1Yr2012_21
  acsByPUMA_C <- DDP.cachedACSa5ByPUMA srcWindow cachedSrc 2021
  let outerKey :: ([GT.StateAbbreviation, GT.PUMA] F.⊆ qs) => F.Record qs -> F.Record [GT.StateAbbreviation, GT.PUMA]
      outerKey = F.rcast
      catKeyO :: (ksO F.⊆ qs) => F.Record qs -> F.Record ksO
      catKeyO = F.rcast
      catKeyM :: (ksM F.⊆ qs) => F.Record qs -> F.Record ksM
      catKeyM = F.rcast
      count = view DT.popCount
      countS = Sum . realToFrac . count
      takeEach n = fmap snd . List.filter ((== 0) . flip mod n . fst) . zip [0..]
      thin = maybe id takeEach thinM
      dataCacheKey = cacheRoot <> "/projModelData.bin"
  let projDataF acsByPUMA = do
        let pdByPUMA = FL.fold (productDistributionFld ms DTP.sumLens outerKey catKeyO countS) acsByPUMA
        projRows <- K.knitEither $ rowsWithProjectedDiffs mc.projVecs pdByPUMA outerKey catKeyO $ thin $ FL.fold FL.list acsByPUMA
        pure $ ProjData (modelNumNullVecs mc) (SBB.rowLength mc.predDMR) projRows
  when clearCaches $ BRCC.clearIfPresentD dataCacheKey
  modelData_C <- BRCC.retrieveOrMakeD (cacheRoot <> "/projModelData.bin") acsByPUMA_C projDataF
  let meanSDFld :: FL.Fold Double (Double, Double) = (,) <$> FL.mean <*> FL.std
      meanSDFlds :: Int -> FL.Fold [Double] [(Double, Double)]
      meanSDFlds m = traverse (\n -> FL.premap (List.!! n) meanSDFld) [0..(m - 1)]
  modelData <- K.ignoreCacheTime modelData_C
  let meanSDs = FL.fold (FL.premap (\(ProjDataRow _ v) -> VS.toList v) $ meanSDFlds (modelNumNullVecs mc)) $ pdRows modelData
  K.logLE K.Info $ "meanSDs=" <> show meanSDs
  states <-  FL.fold (FL.premap (view GT.stateAbbreviation) FL.set) <$> K.ignoreCacheTime acsByPUMA_C
  (dw, code) <-  S.dataWranglerAndCode modelData_C (pure ())
                 (stateGroupBuilder (view GT.stateAbbreviation)  (S.toList states))
                 (const $ pure ())
                 (\projDataRTT _ -> projModel rc catKeyM count predF mc projDataRTT)

  let nNullVecs = modelNumNullVecs mc
      unwraps = (\n -> S.UnwrapExpr ("matrix(ncol="
                                       <> show nNullVecs
                                       <> ", byrow=TRUE, unlist(jsonData $ nvp_ProjectionData))[,"
                                       <> show n <> "]") ("obsNVP_" <> show n))
                <$> [1..nNullVecs]
  res_C <- S.runModel' @BRCC.SerializerC @BRCC.CacheData
           cacheDirE
           (Right runnerInputNames)
           (Just $ S.StanMCParameters 4 4 (Just 1000) (Just 1000) Nothing Nothing (Just 1))
           dw
           code
           S.DoNothing
           (S.ShinyStan unwraps) --(S.Both [S.UnwrapNamed "successes" "yObserved"])
           modelData_C
           (pure ())
  K.logLE K.Info "projModel run complete."
  pure res_C

newtype PModel1 a = PModel1 { pdLogDensity :: a }
  deriving stock (Show, Functor, Foldable, Traversable, Generic)
  deriving anyclass Flat.Flat

instance Applicative PModel1 where
  pure = PModel1
  (PModel1 f) <*> (PModel1 x) = PModel1 (f x)

data PModel0 a = PModel0
  deriving stock (Show, Functor, Foldable, Traversable, Generic)
  deriving anyclass Flat.Flat

instance Applicative PModel0 where
  pure _ = PModel0
  PModel0 <*> PModel0 = PModel0


designMatrixRow1 :: S.DesignMatrixRow (PModel1 Double)
designMatrixRow1 = S.DesignMatrixRow "Model1" [S.DesignMatrixRowPart "logDensity" 1 (VU.singleton . pdLogDensity)]


designMatrixRow0 :: S.DesignMatrixRow (PModel0 Double)
designMatrixRow0 = S.DesignMatrixRow "PModel0" []

designMatrixRow_1 :: S.DesignMatrixRow (F.Record '[DT.Education4C])
designMatrixRow_1 = S.DesignMatrixRow "Base" [cRP]
  where
    cRP = S.DesignMatrixRowPart "Ones" 1 (const $ VU.singleton 1) -- for pure (state-level) alpha
--    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)


designMatrixRow_1_E :: S.DesignMatrixRow (F.Record '[DT.Education4C])
designMatrixRow_1_E = S.DesignMatrixRow "E" [cRP, eRP]
  where
    cRP = S.DesignMatrixRowPart "Ones" 1 (const $ VU.singleton 1) -- for pure (state-level) alpha
    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)


designMatrixRow_1_S_E_R :: S.DesignMatrixRow (F.Record [DT.SexC, DT.Education4C, DT.Race5C])
designMatrixRow_1_S_E_R = S.DesignMatrixRow "S_E_R" [cRP, sRP, eRP, rRP]
  where
    cRP = S.DesignMatrixRowPart "Ones" 1 (const $ VU.singleton 1) -- for pure (state-level) alpha
    sRP = S.boundedEnumRowPart Nothing "Sex" (view DT.sexC)
    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)
    rRP = S.boundedEnumRowPart (Just DT.R5_WhiteNonHispanic) "Race" (view DT.race5C)
