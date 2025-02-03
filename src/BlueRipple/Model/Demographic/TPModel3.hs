{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
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

module BlueRipple.Model.Demographic.TPModel3
  (
    module BlueRipple.Model.Demographic.TPModel3
  )
where

import qualified BlueRipple.Data.CachingCore as BRKU
import qualified BlueRipple.Model.Demographic.EnrichData as DED
import qualified BlueRipple.Model.Demographic.MarginalStructure as DMS
import qualified BlueRipple.Model.Demographic.TableProducts as DTP
import qualified BlueRipple.Model.StanTools as MST

import qualified BlueRipple.Data.CachingCore as BRCC
import qualified BlueRipple.Data.Keyed as BRK
import qualified BlueRipple.Data.Types.Demographic as DT
import qualified BlueRipple.Data.Types.Geographic as GT

import qualified Knit.Report as K hiding (elements)

import qualified Control.MapReduce.Simple as MR

import qualified Control.Foldl as FL
import qualified Data.IntMap.Strict as IM
import qualified Data.Map.Strict as M
import qualified Data.Set as S

import qualified Data.List as List
import qualified Frames as F
import qualified Frames.Melt as F
import qualified Frames.Serialize as FS
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric
import qualified Data.Vinyl as Vinyl
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Control.Lens (Lens', view, over, (^.), _2)
import GHC.TypeLits (Symbol)

import qualified Stan as S
import qualified Stan.BuildingBlocks as SBB (rowLength)
import Stan (TypedList(..))
import Stan.Operators
import qualified CmdStan as CS

import qualified Flat
import Flat.Instances.Vector ()

data ProjDataRow outerK =
  ProjDataRow { pdKey :: outerK, pdCovariates :: VS.Vector Double, pdCoeff :: Double }

-- NB: nullVecs we use are not the ones from SVD but a subset of a rotation of those via
-- the eigenvectors of the covariance
nullVecProjectionsModelDataFld ::  forall outerK k row w .
                                   (Ord outerK)
                               => Lens' w Double
                               -> DMS.MarginalStructure w k
                               -> DTP.NullVectorProjections k
                               -> (row -> outerK)
                               -> (row -> k)
                               -> (row -> w)
                               -> FL.Fold row (VS.Vector Double) -- covariates
                               -> FL.Fold row [(outerK, VS.Vector Double, VS.Vector Double)]
nullVecProjectionsModelDataFld wl ms nvps outerKey catKey datF datFold = case ms of
  DMS.MarginalStructure _ _ -> MR.mapReduceFold
                               MR.noUnpack
                               (MR.assign outerKey id)
                               (MR.foldAndLabel innerFld (\ok (d, v) -> (ok, d, v)))
    where
      projFld = DTP.diffProjectionsFromJointFld ms wl (DTP.fullToProj nvps) catKey datF
      innerFld = (,) <$> datFold <*> projFld

newtype NVProjectionRowData ks =
  NVProjectionRowData { unNVProjectionRowData :: (F.Record ks, VS.Vector Double, VS.Vector Double)} deriving stock (Generic)

instance (FS.RecFlat ok, Vinyl.RMap ok) => Flat.Flat (NVProjectionRowData ok) where
  size (NVProjectionRowData (r, cvs, ps)) = Flat.size (FS.toS r, cvs, ps)
  encode (NVProjectionRowData (r, cvs, ps)) = Flat.encode (FS.toS r, cvs, ps)
  decode = fmap (\(sr, cvsL, psL) -> NVProjectionRowData (FS.fromS sr, cvsL, psL)) Flat.decode

toProjDataRow :: Int -> (k, VS.Vector Double, VS.Vector Double) -> ProjDataRow k
toProjDataRow n (k, covariates, projections) = ProjDataRow k covariates (projections VS.! n)

cwdInnerFld :: (Ord k, BRK.FiniteSet k)
            => (F.Record rs -> k)
            -> (F.Record rs -> DMS.CellWithDensity)
            -> FL.Fold (F.Record rs) [DMS.CellWithDensity]
cwdInnerFld keyF datF = fmap M.elems $ marginalVecFld keyF
  where
    marginalVecFld f = FL.premap (\r -> (f r, datF r)) (DMS.normalizeAndFillMapFld DMS.cwdWgtLens)
-- fmap (VS.fromList . fmap DMS.cwdWgt . M.elems)

bLogit :: Double -> Double -> Double
bLogit eps x
  | x < eps = f eps
  | x > 1 - eps = f (1 - eps)
  | otherwise = f x
  where
    f y = Numeric.log (y / (1 - y))

cwdListToLogitVec :: [DMS.CellWithDensity] -> VS.Vector Double
cwdListToLogitVec = VS.fromList . fmap (bLogit 1e-10 . DMS.cwdWgt)

cwdListToLogPWDensity :: [DMS.CellWithDensity] -> Double
cwdListToLogPWDensity = --FL.fold (safeLogDiv <$> FL.premap (\cw -> DMS.cwdWgt cw * DMS.cwdDensity cw) FL.sum <*> FL.premap DMS.cwdWgt FL.sum)
  posLog . snd . FL.fold (DT.densityAndPopFld' DT.Geometric (const 1) DMS.cwdWgt DMS.cwdDensity)

posLog :: Double -> Double
posLog z = if z < 1 then 0 else Numeric.log z

{-
safeLogDiv :: Double -> Double -> Double
safeLogDiv x y = if y < 1e-10 then 0
                 else let z = x / y
                      in if z < 1 then 0 else Numeric.log z
-}
{-
cwdCovariatesFld :: (Ord k, BRK.FiniteSet k)
                 => (F.Record rs -> k)
                 -> (F.Record rs -> DMS.CellWithDensity)
                 -> FL.Fold (F.Record rs) (VS.Vector Double)
cwdCovariatesFld keyF datF = fmap (\cws -> VS.concat [VS.singleton (safeLog $ cwdListToPWDensity cws), cwdListToLogitVec cws]) $ cwdInnerFld keyF datF
-}

mergeInnerFlds :: [FL.Fold (F.Record rs) (VS.Vector Double)] -> FL.Fold (F.Record rs) (VS.Vector Double)
mergeInnerFlds = fmap VS.concat . sequenceA

dmr :: Text -> Int -> S.DesignMatrixRow (VS.Vector Double)
dmr t n = S.DesignMatrixRow t [S.DesignMatrixRowPart t n VU.convert]

-- NB: nullVecs we use are not the ones from SVD but a subset of a rotation of those via
-- the eigenvectors of the covariance
nullVecProjectionsModelDataFldCheck ::  forall outerK k row w .
                                        (Ord outerK)
                                    => Lens' w Double
                                    -> DMS.MarginalStructure w k
                                    -> DTP.NullVectorProjections k
                                    -> (row -> outerK)
                                    -> (row -> k)
                                    -> (row -> w)
                                    -> FL.Fold row (VS.Vector Double) -- covariates
                                    -> FL.Fold row [(outerK
                                                    , VS.Vector Double -- covariates
                                                    , VS.Vector Double -- nvProjections
                                                    , [(k, w)] -- product table
                                                    , [(k, w)] -- orig table
                                                    )]
nullVecProjectionsModelDataFldCheck wl ms nvps outerKey catKey datF datFold = case ms of
  DMS.MarginalStructure _ _ptFld -> MR.mapReduceFold
                                    MR.noUnpack
                                    (MR.assign outerKey id)
                                    (MR.foldAndLabel innerFld (\ok (d, (v, pKWs, oKWs)) -> (ok, d, v, pKWs, oKWs)))
    where
--      pcF :: [(k, w)] -> [(k, w)]
--      pcF =  FL.fold ptFld
      results kws = let kws' = DMS.normalize (_2 . wl) kws -- normalized original probs
--                        n = FL.fold (FL.premap (view $ _2 . wl) FL.sum) kws
                    in (DTP.diffProjectionsFromJointKeyedList ms wl (DTP.fullToProj nvps) kws'
                      , FL.fold (DMS.marginalProductFromJointFld wl ms) kws --fmap (over (_2 . wl) (* n)) $ pcF kws -- product at original size
                      , kws
                      )
      projFld = fmap results $ DTP.labeledRowsToKeyedListFld catKey datF
      innerFld = (,) <$> datFold <*> projFld

data ProjData outerK =
  ProjData
  {
    pdNPredictrs :: Int
  , pdRows :: [ProjDataRow outerK]
  }

modelIDT :: S.InputDataType S.ModelDataT (ProjData outerK)
modelIDT = S.ModelData

type ProjDataRTT a = S.RowTypeTag (ProjDataRow a)
{-
newtype RecordKey ks = RecordKey (F.Record ks)

instance Text.Show.Show (F.Record ks) => Show (RecordKey ks) where
  show (RecordKey r) = "RecordKey " <> show r

instance Eq (F.Record ks) => Eq (RecordKey ks) where
  RecordKey r1 == RecordKey r2 = r1 == r2

instance Ord (F.Record ks) => Ord (RecordKey ks) where
  compare (RecordKey r1) (RecordKey r2) = compare r1 r2

instance FS.RecFlat ks => Flat.Flat (RecordKey ks) where
  size (RecordKey k) = Flat.size $ FS.toS k
  encode (RecordKey k) = Flat.encode $ FS.toS k
  decode = fmap (RecordKey . FS.fromS) $ Flat.decode

instance BRK.FiniteSet (F.Record k) => BRK.FiniteSet (RecordKey k) where
  elements = RecordKey <$> BRK.elements @(F.Record k)
-}
data ComponentPredictor g =
  ComponentPredictor { mrGeoAlpha :: Map g Double, mrSI :: V.Vector (Double, Double) }
  | ComponentMean { cpMean :: Double }
  deriving stock (Generic)

deriving anyclass instance (Ord g, Flat.Flat g) => Flat.Flat (ComponentPredictor g)

data Predictor k g = Predictor {predPTD :: DTP.ProjectionsToDiff k, predCPs :: [ComponentPredictor g] } deriving stock (Generic)

deriving anyclass instance (Ord g, Ord k, BRK.FiniteSet k, Flat.Flat g) => Flat.Flat (Predictor k g)

{-
mapPredictor :: DMS.IsomorphicKeys a b -> Predictor a g -> Predictor b g
mapPredictor ik@(DMS.IsomorphicKeys abF _) (Predictor nvps predCPs) =
  Predictor (DTP.mapNullVectorProjections ik nvps) (M.mapKeys abF $ fmap (mapCPKey abF) predCPs)
-}

modelResultNVP :: (Show g, Ord g)
               => ComponentPredictor g
               -> g
               -> VS.Vector Double
               -> Either Text Double
modelResultNVP cp g md = case cp of
  ComponentPredictor gam siv -> do
    geoAlpha <- maybeToRight ("geoAlpha lookup failed for gKey=" <> show g <> ". mrGeoAlpha=" <> show gam) $ M.lookup g gam
    let
        applyOne x (b, m) = b * (x - m)
        beta = V.sum $ V.zipWith applyOne (VS.convert md) siv
    pure $ geoAlpha + beta
  ComponentMean x -> pure x


modelResultNVPs :: (Show g, Ord g) => Predictor k g -> g -> VS.Vector Double -> Either Text [Double]
modelResultNVPs p g cvs = traverse (\mr -> modelResultNVP mr g cvs) $ predCPs p

viaNearestOnSimplex :: DTP.ProjectionsToDiff k -> VS.Vector Double -> VS.Vector Double -> K.Sem r (VS.Vector Double)
viaNearestOnSimplex ptd projWs prodV = do
  let n = VS.sum prodV
  pure $ VS.map (* n) $ DTP.projectToSimplex $ DTP.applyNSPWeights ptd projWs (VS.map (/ n) prodV)

-- NB: This function assumes you give it an ordered and complete list of (k, w) pairs
predictedJoint :: forall g k w r . (Show g, Ord g, K.KnitEffects r)
               => DTP.OptimalOnSimplexF r --(DTP.ProjectionsToDiff k -> VS.Vector Double -> VS.Vector Double -> K.Sem r (VS.Vector Double))
               -> Lens' w Double
               -> Predictor k g
               -> g
               -> VS.Vector Double
               -> [(k, w)]
               -> K.Sem r [(k, w)]
predictedJoint onSimplexM wgtLens p gk covariates keyedProduct = do
  let --n = FL.fold (FL.premap (view wgtLens . snd) FL.sum) keyedProduct
      prodV = VS.fromList $ fmap (view wgtLens . snd) keyedProduct

  nvpsPrediction <- K.knitEither $ VS.fromList <$> modelResultNVPs p gk covariates

  onSimplexWgts <- onSimplexM (predPTD p) nvpsPrediction prodV --wgts DTP.projectToSimplex $ DTP.applyNSPWeights (predNVP p) nvpsPrediction (VS.map (/ n) prodV)
--      newWeights = VS.map (* n) onSimplex
  let f (newWgt, (k, w)) = (k, over wgtLens (const newWgt) w)
      predictedTable = fmap f $ zip (VS.toList onSimplexWgts) keyedProduct
      predV = VS.fromList $ fmap (view wgtLens . snd) predictedTable
      checkV = DTP.nvpConstraints (DTP.nullVectorProjections $ predPTD p) LA.#> (predV - prodV)
  K.logLE (K.Debug 1)
    $ "Region=" <> show gk
    <> "\ncovariates = " <> DED.prettyVector covariates
    <> "\npredicted projections = " <> DED.prettyVector nvpsPrediction
    <> "\npredicted projections (onSimplex) = " <> DED.prettyVector onSimplexWgts
    <> "\npredicted result = " <> DED.prettyVector predV
    <> "\nC * (predicted - product) = " <> DED.prettyVector checkV
  pure predictedTable

stateG :: S.GroupTypeTag Text
stateG = S.GroupTypeTag "State"

stateGroupBuilder :: forall f outerK . (Foldable f, Typeable outerK)
                  => (outerK -> Text) -> f Text -> S.StanDataBuilderEff S.ModelDataT (ProjData outerK) (ProjDataRTT outerK)
stateGroupBuilder saF states = do
  projData <- S.addData "ProjectionData" (modelIDT @outerK) (S.ToFoldable pdRows)
  S.addGroupIndexForData (modelIDT @outerK) stateG projData $ S.makeIndexFromFoldable show (saF . pdKey) states
  S.addGroupIntMapForData (modelIDT @outerK) stateG projData $ S.dataToIntMapFromFoldable (saF . pdKey) states
  pure projData

data ProjModelData outerK =
  ProjModelData
  {
    projDataTag :: ProjDataRTT outerK
  , nPredictorsE :: S.IntE
  , predictorsE :: S.MatrixE
  , projectionsE :: S.VectorE
  }

data AlphaModel = AlphaSimple | AlphaHierCentered | AlphaHierNonCentered deriving stock (Show)
data ThetaModel = ThetaSimple | ThetaHierarchical

alphaModelText :: AlphaModel -> Text
alphaModelText AlphaSimple = "AS"
alphaModelText AlphaHierCentered = "AHC"
alphaModelText AlphaHierNonCentered = "AHNC"

thetaModelText :: ThetaModel -> Text
thetaModelText ThetaSimple = "TS"
thetaModelText ThetaHierarchical = "TH"

data Distribution = NormalDist | CauchyDist | StudentTDist

distributionText :: Distribution -> Text
distributionText CauchyDist = "cauchy"
distributionText NormalDist = "normal"
distributionText StudentTDist = "studentT"

data ModelConfig =
  ModelConfig
  {
    standardizeNVs :: Bool
  , designMatrixRow :: S.DesignMatrixRow (VS.Vector Double)
  , alphaModel :: AlphaModel
  , thetaModel :: ThetaModel
  , distribution :: Distribution
  }

data MeanOrModel = Mean | Model !ModelConfig

momText :: MeanOrModel -> Text
momText Mean = "mean"
momText (Model mc) = modelText mc

modelText :: ModelConfig -> Text
modelText (ModelConfig _ dmr' am tm d) = distributionText d <> "_" <> dmr'.dmName <> "_" <> alphaModelText am <> "_" <> thetaModelText tm

dataText :: ModelConfig -> Text
dataText (ModelConfig _ dmr' _ _ _) = dmr'.dmName

projModelData :: forall outerK . Typeable outerK
              => ModelConfig
              -> ProjDataRTT outerK
              -> S.StanModelBuilderEff (ProjData outerK) () (ProjModelData outerK)
projModelData mc projData = do
--  projData <- S.dataSetTag @(ProjDataRow outerK) S.ModelData "ProjectionData"
--  let projMER :: S.MatrixRowFromData (ProjDataRow outerK) --(outerK, md Double, VS.Vector Double)
--      projMER = S.MatrixRowFromData "nvp" Nothing (modelNumNullVecs mc) (\(_, _, v) -> VU.convert v)
  pmE <- S.addRealData (modelIDT @outerK) projData "projection" Nothing Nothing pdCoeff
  let (_, nPredictorsE') = S.designMatrixColDimBinding mc.designMatrixRow Nothing
  dmE <- if SBB.rowLength mc.designMatrixRow > 0
         then S.addDesignMatrix (modelIDT @outerK) projData (contramap pdCovariates mc.designMatrixRow) Nothing
         else pure $ S.namedE "ERROR" S.SMat -- this shouldn't show up in stan code at all
  pure $ ProjModelData projData nPredictorsE' dmE pmE

-- given K null vectors, S states, and D predictors
-- alpha, theta, sigma
-- alpha is a K row-vector or S x K matrix
data Alpha = SimpleAlpha (S.Parameter S.EReal) | HierarchicalAlpha (S.Parameter S.ECVec)
-- theta is a D x K matrix (or Nothing)
newtype Theta = Theta (Maybe (S.Parameter S.ECVec))
-- sigma is a K row-vector
newtype Sigma = Sigma {unSigma :: S.Parameter S.EReal}

newtype Nu = Nu { unNu :: S.Parameter S.EReal }

data ProjModelParameters where
  NormalProjModelParameters :: Alpha -> Theta -> Sigma -> ProjModelParameters
  CauchyProjModelParameters :: Alpha -> Theta -> Sigma -> ProjModelParameters
  StudentTProjModelParameters :: Alpha -> Theta -> Sigma -> Nu -> ProjModelParameters

paramTheta :: ProjModelParameters -> Theta
paramTheta (NormalProjModelParameters _ t _) = t
paramTheta (CauchyProjModelParameters _ t _) = t
paramTheta (StudentTProjModelParameters _ t _ _) = t

projModelParameters :: ModelConfig -> ProjModelData outerK -> S.StanModelBuilderEff (ProjData outerK) () ProjModelParameters
projModelParameters mc pmd = do
  let stdNormalDWA :: (S.TypeOneOf t [S.EReal, S.ECVec, S.ERVec], S.GenSType t) => S.DensityWithArgs t
      stdNormalDWA = S.DensityWithArgs S.std_normal TNil
      numPredictors = SBB.rowLength mc.designMatrixRow

  sigma <-  Sigma
             <$> S.simpleParameterWA
             (S.NamedDeclSpec "sigma" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)
             stdNormalDWA
  let nStatesE = S.groupSizeE stateG
      hierAlphaNDS = S.NamedDeclSpec "alpha" $ S.vectorSpec nStatesE
      hierAlphaPs = do
        muAlphaP <- S.simpleParameterWA
                    (S.NamedDeclSpec "muAlpha" S.realSpec)
                    stdNormalDWA
        sigmaAlphaP <-  S.simpleParameterWA
                        (S.NamedDeclSpec "sigmaAlpha" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)
                        stdNormalDWA
        pure (muAlphaP :> sigmaAlphaP :> TNil)
  alpha <- case mc.alphaModel of
    AlphaSimple -> do
      fmap SimpleAlpha
        $ S.simpleParameterWA
        (S.NamedDeclSpec "alpha" S.realSpec)
        stdNormalDWA
    AlphaHierCentered -> do
      alphaPs <- hierAlphaPs
      fmap HierarchicalAlpha
        $ S.addBuildParameter
        $ S.UntransformedP hierAlphaNDS [] alphaPs
        $ \(muAlphaE :> sigmaAlphaE :> TNil) m
          -> S.addStmt $ S.sample m S.normalS (muAlphaE :> sigmaAlphaE :> TNil)
    AlphaHierNonCentered -> do
      alphaPs <- hierAlphaPs
      let rawNDS = S.NamedDeclSpec (S.declName hierAlphaNDS <> "_raw") $ S.decl hierAlphaNDS
      rawAlphaP <- S.simpleParameterWA rawNDS stdNormalDWA
      fmap HierarchicalAlpha
        $ S.addBuildParameter
        $ S.TransformedP hierAlphaNDS []
        (rawAlphaP :> alphaPs) S.TransformedParametersBlock
        (\(rawE :> muAlphaE :> muSigmaE :> TNil) -> S.DeclRHS $ muAlphaE `S.plusE` (muSigmaE `S.timesE` rawE))
        TNil (\_ _ -> pure ())
  -- for now all the thetas are iid normal, but perchance hierarchical
  theta <- if numPredictors == 0
           then pure (Theta Nothing)
           else case mc.thetaModel of
                  ThetaSimple ->
                    (Theta . Just)
                    <$> S.simpleParameterWA
                    (S.NamedDeclSpec "theta" $ S.vectorSpec pmd.nPredictorsE)
                    stdNormalDWA
                  ThetaHierarchical -> do
                    sigmaThetaP <-  S.simpleParameterWA
                      (S.NamedDeclSpec "sigmaTheta" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)
                      stdNormalDWA
                    fmap (Theta . Just)
                      $ S.addBuildParameter
                      $ S.UntransformedP
                      (S.NamedDeclSpec "theta" $ S.vectorSpec pmd.nPredictorsE) [] (sigmaThetaP :> TNil)
                      $ \(sigmaThetaE :> TNil) m -> S.addStmt $ S.sample m S.normalS (S.realE 0 :> sigmaThetaE :> TNil)
  case mc.distribution of
    NormalDist -> pure $ NormalProjModelParameters alpha theta sigma
    CauchyDist -> pure $ CauchyProjModelParameters alpha theta sigma
    StudentTDist -> do
--      let kVectorOf x = S.functionE S.rep_row_vector (S.realE x :> pmd.nNullVecsE :> TNil)
      nu <-  fmap Nu
             $ S.simpleParameterWA
             (S.NamedDeclSpec "nu" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)
             (S.DensityWithArgs S.gamma (S.realE 2 :> S.realE 0.1 :> TNil))
      pure $ StudentTProjModelParameters alpha theta sigma nu

data RunConfig = RunConfig { nvIndex :: Int, rcIncludePPCheck :: Bool, rcIncludeLL :: Bool, statesM :: Maybe (Text, [Text]) }

-- not returning anything for now
projModel :: (Typeable outerK)
          => RunConfig
          -> ModelConfig
          -> ProjDataRTT outerK
          -> S.StanModelBuilderEff (ProjData outerK) () ()
projModel rc mc projData = do
  mData <- projModelData mc projData
  mParams <- projModelParameters mc mData
  let betaNDS = S.NamedDeclSpec "beta" $ S.vectorSpec mData.nPredictorsE
      nRowsE = S.dataSetSizeE mData.projDataTag
--      fstI x k = S.sliceE S.s0 k x
      pExpr = S.parameterExpr

--      sndI x k = S.sliceE S.s1 k x
--      loopNVs = S.for "k" (S.SpecificNumbered (S.intE 1) mData.nNullVecsE)
  (predM, _centerF, _mBeta) <- case paramTheta mParams of
    Theta (Just thetaP) -> do
      (centeredPredictorsE, centerF) <- S.centerDataMatrix S.DMCenterOnly mData.predictorsE Nothing "DM"
      (dmQ, _, _, mBeta) <- S.thinQR centeredPredictorsE "DM" $ Just (pExpr thetaP, betaNDS)
      pure (dmQ, centerF, mBeta)
    Theta Nothing -> pure (S.namedE "ERROR" S.SMat, \_ x _ -> pure x, Nothing)
  (stdNVP, inverseF) <- case mc.standardizeNVs of
    True -> do
      sdP <- S.addBuildParameter
             $ S.TransformedDataP
             $ S.TData
             (S.NamedDeclSpec "nvpSD" S.realSpec) []
             TNil (\_ -> S.DeclRHS $ S.sd mData.projectionsE)
      stdNVP <- S.inBlock S.SBTransformedData $ S.addFromCodeWriter
                $ S.declareRHSNW (S.NamedDeclSpec "stdNVP" $ S.vectorSpec nRowsE)
                $ mData.projectionsE |/| pExpr sdP
      let inverse :: (t ~ S.BinaryResultT S.BMultiply S.EReal t) => S.UExpr t -> S.UExpr t
          inverse psCol = pExpr sdP `S.timesE` psCol
      pure (stdNVP, inverse)
    False -> pure (mData.projectionsE, id)

  -- model
  let reIndexByState = S.indexE S.s0 (S.dataByGroupIndexE mData.projDataTag stateG)
      muE :: Alpha -> Theta -> S.VectorE
      muE a t =  case a of
       SimpleAlpha alphaP -> case t of
         Theta Nothing -> S.rep_vector (pExpr alphaP) nRowsE
         Theta (Just thetaP) -> pExpr alphaP |+| (predM |*| pExpr thetaP)
       HierarchicalAlpha alpha -> case t of
         Theta Nothing -> reIndexByState $ pExpr alpha
         Theta (Just thetaP) -> reIndexByState (pExpr alpha) |+| (predM |*| pExpr thetaP)
      sigmaE :: Sigma -> S.VectorE
      sigmaE s = S.rep_vector (pExpr (unSigma s)) nRowsE

  let ppF :: ((S.IntE -> S.ExprList xs) -> S.IntE -> S.UExpr S.EReal)
          -> S.CodeWriter (S.IntE -> S.ExprList xs)
          -> S.StanModelBuilderEff (ProjData outerK) () (S.ArrayE S.EReal)
      ppF rngF rngPSCW = S.generatePosteriorPrediction'
                         mData.projDataTag
                         (S.NamedDeclSpec "predProj" $ S.array1Spec nRowsE S.realSpec)
                         rngF
                         rngPSCW
                         (\_ p -> inverseF p)
      llF :: S.StanDist t pts rts
          -> S.CodeWriter (S.IntE -> S.ExprList pts)
          -> S.CodeWriter (S.IntE -> S.UExpr t)
          -> S.StanModelBuilderEff md gq ()
      llF = S.generateLogLikelihood mData.projDataTag

  let (sampleStmtF, pp, ll) = case mParams of
        NormalProjModelParameters a t s ->
          let ssF e = S.sample e S.normal (muE a t :> sigmaE s  :> TNil)
              rF f nE = S.functionE S.normal_rngF (f nE)
              rpF = pure $ \nE -> muE a t !! nE :> sigmaE s !! nE :> TNil
              ll' = llF S.normalDist rpF (pure $ \nE -> stdNVP !! nE)
          in (ssF, ppF rF rpF, ll')
        CauchyProjModelParameters a t s ->
          let ssF e = S.sample e S.cauchy (muE a t :> sigmaE s :> TNil)
              rF f nE = S.functionE S.cauchy_rngF (f nE)
              rpF = pure $ \nE -> muE a t !! nE :> sigmaE s !! nE :> TNil
              ll' = llF S.cauchyDist rpF (pure $ \nE -> stdNVP !! nE)
          in (ssF, ppF rF rpF, ll')
        StudentTProjModelParameters a t s n ->
          let nu :: Nu -> S.VectorE
              nu n' = S.rep_vector (pExpr (unNu n')) nRowsE
              ssF e = S.sample e S.student_t (nu n  :> muE a t :> sigmaE s :> TNil)
              rF f nE = S.functionE S.student_t_rngF (f nE)
              rpF = pure $ \nE -> nu n !! nE :> muE a t !! nE :> sigmaE s !! nE :>  TNil
              ll' = llF S.studentTDist rpF (pure $ \nE -> stdNVP !! nE)
          in (ssF, ppF rF rpF, ll')

  S.inBlock S.SBModel $ S.addFromCodeWriter $ S.addStmt $ sampleStmtF stdNVP
  -- generated quantities
  when rc.rcIncludePPCheck $ void pp
  when rc.rcIncludeLL ll
  pure ()

cwdF :: (F.ElemOf rs DT.PopCount, F.ElemOf rs DT.PWPopPerSqMile) => F.Record rs -> DMS.CellWithDensity
cwdF r = DMS.CellWithDensity (realToFrac $ r ^. DT.popCount) (r ^. DT.pWPopPerSqMile)

model3A5CacheDir :: Text
model3A5CacheDir = "model/demographic/nullVecProjModel3_A5/"

buildProjModelNVPData ::  forall (ks :: [(Symbol, Type)]) rs r .
                          (K.KnitEffects r
                          , BRKU.CacheEffects r
                          , ks F.⊆ rs
                          , F.ElemOf rs GT.PUMA
                          , F.ElemOf rs GT.StateAbbreviation
                          , F.ElemOf rs DT.PopCount
                          , F.ElemOf rs DT.PWPopPerSqMile
                          )
                      => Either Text Text
                      -> Text
                      -> K.ActionWithCacheTime r (F.FrameRec rs)
                      -> K.ActionWithCacheTime r (DTP.NullVectorProjections (F.Record ks))
                      -> DMS.MarginalStructure DMS.CellWithDensity (F.Record ks)
                      -> FL.Fold (F.Record rs) (VS.Vector Double)
                      -> K.Sem r (K.ActionWithCacheTime r [NVProjectionRowData [GT.StateAbbreviation, GT.PUMA]])
buildProjModelNVPData cacheDirE modelId acs_C nvps_C ms datFld = K.wrapPrefix "TPModel3.buildProjModelNVPData" $ do
  nvpDataCacheKey <- BRKU.cacheFromDirE cacheDirE (modelId <> "_nvpRowData.bin")
  let outerKey = F.rcast @[GT.StateAbbreviation, GT.PUMA]
      catKey = F.rcast @ks
      nvpDeps = (,) <$> acs_C <*>  nvps_C
      rawRows acs nvps = NVProjectionRowData <$> FL.fold (nullVecProjectionsModelDataFld DMS.cwdWgtLens ms nvps outerKey catKey cwdF datFld) acs
      logRebuild = K.logLE K.Info "(Re)building data for models of each projection."
  BRKU.retrieveOrMakeD nvpDataCacheKey nvpDeps $ \(x, y) ->  logRebuild >> (pure $ rawRows x y)


runProjModel :: forall (ks :: [(Symbol, Type)]) rs r .
                (K.KnitEffects r
                , BRKU.CacheEffects r
                , ks F.⊆ rs
                , F.ElemOf rs GT.PUMA
                , F.ElemOf rs GT.StateAbbreviation
                , F.ElemOf rs DT.PopCount
                , F.ElemOf rs DT.PWPopPerSqMile
                )
             => Either Text Text
             -> RunConfig
             -> MeanOrModel
             -> K.ActionWithCacheTime r [NVProjectionRowData [GT.StateAbbreviation, GT.PUMA]]
             -> Set Text
             -> K.Sem r (K.ActionWithCacheTime r (ComponentPredictor Text))
runProjModel cacheDirE rc mom projData_C states {- nvps_C ms datFld -} = K.wrapPrefix "TPModel3.runProjModel" $ do
  let statesFilter = maybe id (\(_, sts) -> filter ((`elem` sts) . view GT.stateAbbreviation . pdKey)) rc.statesM
      rowData_C = fmap (statesFilter . fmap (toProjDataRow rc.nvIndex)) $ fmap (fmap unNVProjectionRowData) projData_C
  case mom of
    Mean -> do
      let meanSDFld :: FL.Fold Double (Double, Double) = (,) <$> FL.mean <*> FL.std
      rowData <- K.ignoreCacheTime rowData_C
      let (mean, sd) = FL.fold (FL.premap pdCoeff meanSDFld) $ rowData
      K.logLE K.Info $ "mean=" <> show mean <> "; sd=" <> show sd
      pure $ pure $ ComponentMean mean
    Model mc -> do
      let modelData_C = ProjData (SBB.rowLength mc.designMatrixRow) <$> rowData_C
          dataName = "projectionData_" <> dataText mc <> "_N" <> show rc.nvIndex <> maybe "" fst rc.statesM
      stanDir <- K.liftKnit MST.stanDir >>= K.knitMaybe "runModel: empty stanDir!" . BRCC.insureFinalSlash
      let runnerInputNames = S.RunnerInputNames
                             (stanDir <> "demographic/nullVecProj_M3_A5")
                             (modelText mc)
                             (Just $ S.GQNames "pp" dataName) -- posterior prediction vars to wrap
                             dataName
      (dw, code) <-  S.dataWranglerAndCode modelData_C (pure ())
                     (stateGroupBuilder (view GT.stateAbbreviation)  (S.toList states))
                     (const $ pure ())
                     (\projDataRTT _  -> projModel rc mc projDataRTT)

      let unwraps = [S.UnwrapNamed "projection" "yProjection"]

      res_C <- S.runModel' @BRKU.SerializerC @BRKU.CacheData
               cacheDirE
               (Right runnerInputNames)
               Nothing
               dw
               code
               (projModelResultAction mc) --S.DoNothing -- (stateModelResultAction mcWithId dmr)
               (S.Both unwraps) --(S.Both [S.UnwrapNamed "successes" "yObserved"])
               modelData_C
               (pure ())
      K.logLE K.Info "projModel run complete."
      pure res_C

--NB: parsed summary data has stan indexing, i.e., Arrays start at 1.
projModelResultAction :: forall outerK r .
                         (K.KnitEffects r
                         , Typeable outerK
                         )
                      => ModelConfig
                      -> S.ResultAction (ProjData outerK) () S.DataSetGroupIntMaps S.DataSetGroupIntMaps r () (ComponentPredictor Text)
projModelResultAction mc = S.UseSummary f where
  f summary _ modelDataAndIndexes_C _ = do
    (modelData, resultIndexesE) <- K.ignoreCacheTime modelDataAndIndexes_C
    -- compute means of predictors because model was zero-centered in them
    let nPredictors = SBB.rowLength mc.designMatrixRow
        mdMeansFld = FL.premap (VS.toList . pdCovariates)
                    $ traverse (\n -> FL.premap (List.!! n) FL.mean) [0..(nPredictors - 1)]
        nvpSDFld = FL.premap pdCoeff FL.std
        (mdMeansL, nvpSD) = FL.fold ((,) <$> mdMeansFld <*> nvpSDFld) $ pdRows modelData
        rescaleAlphaBeta x = if mc.standardizeNVs then x * nvpSD else x
    stateIM <- K.knitEither
      $ resultIndexesE >>= S.getGroupIndex (S.RowTypeTag @(ProjDataRow outerK) "ProjectionData") stateG
    let allStates = IM.elems stateIM
        getScalar n = K.knitEither $ S.getScalar . fmap CS.mean <$> S.parseScalar n (CS.paramStats summary)
        getVector n = K.knitEither $ S.getVector . fmap CS.mean <$> S.parse1D n (CS.paramStats summary)
--        getMatrix n = K.knitEither $ fmap CS.mean <$> S.parse2D n (CS.paramStats summary)
    geoMap <- case mc.alphaModel of
      AlphaSimple -> do
        alpha <- getScalar "alpha" -- states by nNullvecs
        pure $ M.fromList $ fmap (, rescaleAlphaBeta alpha) allStates
      _ -> do
        alphaV <- getVector "alpha"
--        let mRowToList row = S.getIndexed alphaV row
        pure $ M.fromList $ fmap (\(stateIdx, stateAbbr) -> (stateAbbr, rescaleAlphaBeta (alphaV V.! (stateIdx - 1)))) $ IM.toList stateIM
    betaSI <- case nPredictors of
      0 -> pure V.empty
      _p -> do
        betaV <- getVector "beta"
        pure $ V.fromList $ zip (fmap rescaleAlphaBeta $ V.toList betaV) mdMeansL
    pure $ ComponentPredictor geoMap betaSI

{-
data ASERModelP a = ASERModelP { mASER_PWLogDensity :: a, mASER_FracOver45 :: a, mASER_FracGrad :: a, mASER_FracOfColor :: a , mASER_FracWNG :: a  }
  deriving stock (Show, Generic)
  deriving anyclass Flat.Flat

aserModelDatFld :: (F.ElemOf rs DT.PWPopPerSqMile
                   , F.ElemOf rs DT.Age4C
                   , F.ElemOf rs DT.Education4C
                   , F.ElemOf rs DT.Race5C
                   , F.ElemOf rs DT.PopCount
                   )
             => FL.Fold (F.Record rs) (ASERModelP Double)
aserModelDatFld = ASERModelP <$> dFld <*> aFld <*> gFld <*> rFld <*> wngFld
  where
    nPeople = realToFrac . view DT.popCount
    dens r = let x = view DT.pWPopPerSqMile r in if x > 1 then Numeric.log x else 0
    wgtFld = FL.premap nPeople FL.sum
    wgtdFld f = safeDiv <$> FL.premap (\r -> nPeople r * f r) FL.sum <*> wgtFld
    dFld = wgtdFld dens
    over45 = (`elem` [DT.A4_45To64, DT.A4_65AndOver]) . view DT.age4C
    grad = (== DT.E4_CollegeGrad) . view DT.education4C
    ofColor = (/= DT.R5_WhiteNonHispanic) . view DT.race5C
    fracFld f = safeDiv <$> FL.prefilter f wgtFld <*> wgtFld
    wng x = not (ofColor x) && not (grad x)
    aFld = fracFld over45
    gFld = fracFld grad
    rFld = fracFld ofColor
    wngFld = fracFld wng

aserModelFuncs :: ModelDataFuncs ASERModelP a
aserModelFuncs = ModelDataFuncs aserModelToList aserModelFromList where
  aserModelToList :: ASERModelP a -> [a]
  aserModelToList (ASERModelP x y z a b) = [x, y, z, a, b]

  aserModelFromList :: [a] -> Either Text (ASERModelP a)
  aserModelFromList as = case as of
    [x, y, z, a, b] -> Right $ ASERModelP x y z a b
    _ -> Left "aserModelFromList: wrong size list given (n /= 5)"

designMatrixRowASER :: S.DesignMatrixRow (ASERModelP Double)
designMatrixRowASER = S.DesignMatrixRow "ASER"
                      [S.DesignMatrixRowPart "logDensity" 1 (VU.singleton . mASER_PWLogDensity)
                      , S.DesignMatrixRowPart "fracOver45" 1 (VU.singleton . mASER_FracOver45)
                      , S.DesignMatrixRowPart "fracGrad" 1 (VU.singleton . mASER_FracGrad)
                      , S.DesignMatrixRowPart "fracOC" 1 (VU.singleton . mASER_FracOfColor)
                      , S.DesignMatrixRowPart "fracWNG" 1 (VU.singleton . mASER_FracWNG)
                      ]


---

data Model1P a = Model1P { m1pPWLogDensity :: a, m1pFracGrad :: a, m1pFracOfColor :: a }
  deriving stock (Show, Generic)
  deriving anyclass Flat.Flat

model1DatFld :: (F.ElemOf rs DT.PWPopPerSqMile
                , F.ElemOf rs DT.Education4C
                , F.ElemOf rs DT.Race5C
                , F.ElemOf rs DT.PopCount
                )
             => FL.Fold (F.Record rs) (Model1P Double)
model1DatFld = Model1P <$> dFld <*> gFld <*> rFld
  where
    nPeople = realToFrac . view DT.popCount
    dens = Numeric.log . view DT.pWPopPerSqMile
    wgtFld = FL.premap nPeople FL.sum
    wgtdFld f = (/) <$> FL.premap (\r -> nPeople r * f r) FL.sum <*> wgtFld
    dFld = wgtdFld dens
    fracFld f = (/) <$> FL.prefilter f wgtFld <*> wgtFld
    gFld = fracFld ((== DT.E4_CollegeGrad) . view DT.education4C)
    rFld = fracFld ((/= DT.R5_WhiteNonHispanic) . view DT.race5C)

model1Funcs :: ModelDataFuncs Model1P a
model1Funcs = ModelDataFuncs model1ToList model1FromList where
  model1ToList :: Model1P a -> [a]
  model1ToList (Model1P x y z) = [x, y, z]

  model1FromList :: [a] -> Either Text (Model1P a)
  model1FromList as = case as of
    [x, y, z] -> Right $ Model1P x y z
    _ -> Left "model1FromList: wrong size list given (n /= 3)"

emptyDM :: S.DesignMatrixRow (Model1P Double)
emptyDM = S.DesignMatrixRow "EDM" []

designMatrixRow1 :: S.DesignMatrixRow (Model1P Double)
designMatrixRow1 = S.DesignMatrixRow "PM1"
                   [S.DesignMatrixRowPart "logDensity" 1 (VU.singleton . m1pPWLogDensity)
                   , S.DesignMatrixRowPart "fracGrad" 1 (VU.singleton . m1pFracGrad)
                   , S.DesignMatrixRowPart "fracOC" 1 (VU.singleton . m1pFracOfColor)
                   ]

data Model2P a = Model2P { m2pPWLogDensity :: a, m2pFracCit :: a, m2pFracGrad :: a, m2pFracOfColor :: a }
  deriving stock (Show, Generic)
  deriving anyclass Flat.Flat

safeDiv :: Double -> Double -> Double
safeDiv x y = if y /= 0 then x / y else 0
{-# INLINE safeDiv #-}

model2DatFld :: (F.ElemOf rs DT.PWPopPerSqMile
                , F.ElemOf rs DT.CitizenC
                , F.ElemOf rs DT.Education4C
                , F.ElemOf rs DT.Race5C
                , F.ElemOf rs DT.PopCount
                )
             => FL.Fold (F.Record rs) (Model2P Double)
model2DatFld = Model2P <$> dFld <*> cFld <*> gFld <*> rFld
  where
    nPeople = realToFrac . view DT.popCount
    dens r = let pwd = view DT.pWPopPerSqMile r in if pwd > 1 then Numeric.log pwd else 0
    wgtFld = FL.premap nPeople FL.sum
    wgtdFld f = safeDiv <$> FL.premap (\r -> nPeople r * f r) FL.sum <*> wgtFld
    dFld = wgtdFld dens
    fracFld f = (/) <$> FL.prefilter f wgtFld <*> wgtFld
    cFld = fracFld ((== DT.Citizen) . view DT.citizenC)
    gFld = fracFld ((== DT.E4_CollegeGrad) . view DT.education4C)
    rFld = fracFld ((/= DT.R5_WhiteNonHispanic) . view DT.race5C)

model2Funcs :: ModelDataFuncs Model2P a
model2Funcs = ModelDataFuncs model2ToList model2FromList where
  model2ToList :: Model2P a -> [a]
  model2ToList (Model2P x y z a) = [x, y, z, a]

  model2FromList :: [a] -> Either Text (Model2P a)
  model2FromList as = case as of
    [x, y, z, a] -> Right $ Model2P x y z a
    _ -> Left "model1FromList: wrong size list given (n /= 3)"

designMatrixRow2 :: S.DesignMatrixRow (Model2P Double)
designMatrixRow2 = S.DesignMatrixRow "PM2"
                   [S.DesignMatrixRowPart "logDensity" 1 (VU.singleton . m2pPWLogDensity)
                   , S.DesignMatrixRowPart "fracCit" 1 (VU.singleton . m2pFracCit)
                   , S.DesignMatrixRowPart "fracGrad" 1 (VU.singleton . m2pFracGrad)
                   , S.DesignMatrixRowPart "fracOC" 1 (VU.singleton . m2pFracOfColor)
                   ]
-}
