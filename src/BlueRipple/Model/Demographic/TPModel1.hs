{-# LANGUAGE AllowAmbiguousTypes #-}
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

module BlueRipple.Model.Demographic.TPModel1
  (
    module BlueRipple.Model.Demographic.TPModel1
  )
where

--import qualified BlueRipple.Configuration as BR
import qualified BlueRipple.Data.CachingCore as BRCC
import qualified BlueRipple.Model.Demographic.DataPrep as DDP
import qualified BlueRipple.Model.Demographic.MarginalStructure as DMS
import qualified BlueRipple.Model.Demographic.TableProducts as DTP
import qualified BlueRipple.Model.StanTools as MST



import qualified BlueRipple.Data.Types.Demographic as DT
import qualified BlueRipple.Data.Types.Geographic as GT
import qualified BlueRipple.Data.ACS_PUMS as ACS

import qualified Knit.Report as K

import qualified Control.MapReduce.Simple as MR

import qualified Control.Foldl as FL
import qualified Data.IntMap.Strict as IM
import qualified Data.Map.Strict as M
import qualified Data.Set as Set

import qualified Data.List as List
import qualified Frames as F
import qualified Frames.Melt as F
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU

import Control.Lens (view, _2)
import GHC.TypeLits (Symbol)

import qualified Stan as S
import qualified Stan.BuildingBlocks as SBB (rowLength)
import Stan (TypedList(..))
import Stan.Operators
import qualified CmdStan as CS
{-
import qualified Stan.ModelBuilder as SMB
import qualified Stan.ModelRunner as SMR
import qualified Stan.ModelConfig as SC
import qualified Stan.Parameters as SP
import qualified Stan.RScriptBuilder as SR
import qualified Stan.ModelBuilder.BuildingBlocks as SBB
import qualified Stan.ModelBuilder.DesignMatrix as DM
import qualified Stan.ModelBuilder.TypedExpressions.Types as TE
import qualified Stan.ModelBuilder.TypedExpressions.Statements as TE
import qualified Stan.ModelBuilder.TypedExpressions.Indexing as TEI
import qualified Stan.ModelBuilder.TypedExpressions.Operations as TEO
import qualified Stan.ModelBuilder.TypedExpressions.DAG as DAG
import qualified Stan.ModelBuilder.TypedExpressions.StanFunctions as SF
import Stan.ModelBuilder.TypedExpressions.TypedList (TypedList(..))
-}
import qualified Flat

-- NB: nullVecs we use are not the ones from SVD but a subset of a rotation of those via
-- the eigenvectors of the covariance
nullVecProjectionsModelDataFld ::  forall outerK k row md .
                                   (Ord outerK)
                               => DMS.MarginalStructure (Sum Double) k
                               -> DTP.NullVectorProjections k
                               -> (row -> outerK)
                               -> (row -> k)
                               -> (row -> LA.R)
                               -> FL.Fold row md
                               -> FL.Fold row [(outerK, md, LA.Vector LA.R)]
nullVecProjectionsModelDataFld ms nvps outerKey catKey count datFold = case ms of
  DMS.MarginalStructure _ _ -> MR.mapReduceFold
                                   MR.noUnpack
                                   (MR.assign outerKey id)
                                   (MR.foldAndLabel innerFld (\ok (d, v) -> (ok, d, v)))
    where
      projFld = DTP.diffProjectionsFromJointFld ms DTP.sumLens (DTP.fullToProj nvps) catKey (Sum . count)
      innerFld = (,) <$> datFold <*> projFld

-- NB: nullVecs we use are not the ones from SVD but a subset of a rotation of those via
-- the eigenvectors of the covariance
nullVecProjectionsModelDataFldCheck ::  forall outerK k row md .
                                        (Ord outerK)
                                    => DMS.MarginalStructure (Sum Double) k
                                    -> DTP.NullVectorProjections k
                                    -> (row -> outerK)
                                    -> (row -> k)
                                    -> (row -> Double)
                                    -> FL.Fold row md
                                    -> FL.Fold row [(outerK, md, LA.Vector Double, LA.Vector Double, LA.Vector Double)]
nullVecProjectionsModelDataFldCheck ms nvps outerKey catKey count datFold = case ms of
  DMS.MarginalStructure _ ptFld -> MR.mapReduceFold
                                   MR.noUnpack
                                   (MR.assign outerKey id)
                                   (MR.foldAndLabel innerFld (\ok (d, (v, pv, nv)) -> (ok, d, v, pv, nv)))
    where
--      allKs :: Set k = BRK.elements
      pcF :: [(k, Sum Double)] -> VS.Vector Double
      pcF =  VS.fromList . fmap (getSum . snd) . FL.fold ptFld
      results kws = let kws' = DMS.normalize (_2 . DTP.sumLens) kws --DTP.normalizedVec v -- normalized original probs
                        nSum = FL.fold (FL.premap (view $ _2 . DTP.sumLens) FL.sum) kws --VS.sum v -- num of people
                    in (DTP.diffProjectionsFromJointKeyedList ms DTP.sumLens (DTP.fullToProj nvps) kws' -- DTP.fullToProj nvps (w - pcF w)
                       , VS.map (* nSum) (pcF kws')
                       , VS.fromList $ fmap (view $ _2 . DTP.sumLens) kws
                       )
      projFld = fmap results $ DTP.labeledRowsToKeyedListFld catKey (Sum . count)
      innerFld = (,) <$> datFold <*> projFld


type ProjDataRow outerK md = (outerK, md Double, LA.Vector LA.R)

data ProjData outerK md =
  ProjData
  {
    pdNNullVecs :: Int
  , pdNPredictors :: Int
  , pdRows :: [ProjDataRow outerK md]
  }

modelIDT :: forall outerK md . S.InputDataType S.ModelDataT (ProjData outerK md)
modelIDT = S.ModelData

type ProjDataRTT outerK md = S.RowTypeTag (ProjDataRow outerK md)

data ModelDataFuncs md a = ModelDataFuncs { mdfToList :: md a -> [a], mdfFromList :: [[(a, a)]] -> Either Text (md [(a, a)])}


data ModelResult g (b :: Type -> Type) = ModelResult { mrGeoAlpha :: Map g [Double], mrSI :: b [(Double, Double)] }
  deriving stock (Generic)

deriving stock instance (Show g, Show (b [(Double, Double)])) => Show (ModelResult g b)
deriving anyclass instance (Ord g, Flat.Flat g, Flat.Flat (b [(Double, Double)])) => Flat.Flat (ModelResult g b)

modelResultNVPs :: (Show g, Ord g) => (forall a . ModelDataFuncs md a) -> ModelResult g md -> g -> md Double -> Either Text [Double]
modelResultNVPs mdf mr g md = do
    geoAlphas <- maybeToRight ("geoAlpha lookup failed for gKey=" <> show g) $ M.lookup g mr.mrGeoAlpha
    let mdL = mdf.mdfToList md
        mrBetaL = mdf.mdfToList mr.mrSI
        applyOne x (b, m) = b * (x - m)
        applyToList x = fmap (applyOne x)
        eachBetaL = zipWith applyToList mdL mrBetaL
        betaL = fmap (FL.fold FL.sum) $ transp eachBetaL
    pure $ zipWith (+) geoAlphas betaL

stateG :: S.GroupTypeTag Text
stateG = S.GroupTypeTag "State"

stateGroupBuilder :: forall f outerK md . (Foldable f, Typeable outerK, Typeable md)
                  => (outerK -> Text) -> f Text -> S.StanDataBuilderEff S.ModelDataT (ProjData outerK md) (ProjDataRTT outerK md)
stateGroupBuilder saF states = do
  let ok (x, _, _) = x
  projData <- S.addData "ProjectionData" (modelIDT @outerK @md) (S.ToFoldable pdRows)
  S.addGroupIndexForData (modelIDT @outerK @md) stateG projData $ S.makeIndexFromFoldable show (saF . ok) states
  S.addGroupIntMapForData (modelIDT @outerK @md) stateG projData $ S.dataToIntMapFromFoldable (saF . ok) states
  pure projData

data ProjModelData outerK md =
  ProjModelData
  {
    projDataTag :: ProjDataRTT outerK md
  , nNullVecsE :: S.IntE
  , nPredictorsE :: S.IntE
  , predictorsE :: S.MatrixE
  , projectionsE :: S.MatrixE
  }

data AlphaModel = AlphaSimple | AlphaHierCentered | AlphaHierNonCentered deriving stock (Show)

alphaModelText :: AlphaModel -> Text
alphaModelText AlphaSimple = "AS"
alphaModelText AlphaHierCentered = "AHC"
alphaModelText AlphaHierNonCentered = "AHNC"

data Distribution = NormalDist | CauchyDist | StudentTDist

distributionText :: Distribution -> Text
distributionText CauchyDist = "cauchy"
distributionText NormalDist = "normal"
distributionText StudentTDist = "studentT"

data ModelConfig k md =
  ModelConfig
  {
    projVecs :: DTP.NullVectorProjections k
  , standardizeNVs :: Bool
  , designMatrixRow :: S.DesignMatrixRow (md Double)
  , alphaModel :: AlphaModel
  , distribution :: Distribution
  , mdFuncs :: ModelDataFuncs md Double
  }

modelNumNullVecs :: ModelConfig k md -> Int
modelNumNullVecs mc = fst $ LA.size $ DTP.nvpProj mc.projVecs

modelText :: ModelConfig k md -> Text
modelText mc = distributionText mc.distribution <> "_" <> mc.designMatrixRow.dmName <> "_" <> alphaModelText mc.alphaModel

dataText :: ModelConfig k md -> Text
dataText mc = mc.designMatrixRow.dmName <> "_NV" <> show (modelNumNullVecs mc)

projModelData :: forall md outerK k . (Typeable outerK, Typeable md)
              => ModelConfig k md
              -> ProjDataRTT outerK md
              -> S.StanModelBuilderEff (ProjData outerK md) () (ProjModelData outerK md)
projModelData mc projData = do
--  projData <- S.dataSetTag @(ProjDataRow outerK md) S.ModelData "ProjectionData"
  let projMER :: S.MatrixRowFromData (ProjDataRow outerK md) --(outerK, md Double, VS.Vector Double)
      projMER = S.MatrixRowFromData "nvp" Nothing (modelNumNullVecs mc) (\(_, _, v) -> VU.convert v)
  (pmE, nNullVecsE') <- S.add2dMatrixData (modelIDT @outerK @md) projData projMER Nothing Nothing
--  let nNullVecsE' = S.mrfdColumnsE projMER
  let (_, nPredictorsE') = S.designMatrixColDimBinding mc.designMatrixRow Nothing
  dmE <- if SBB.rowLength mc.designMatrixRow > 0
         then S.addDesignMatrix (modelIDT @outerK @md) projData (contramap (\(_, md, _) -> md) mc.designMatrixRow) Nothing
         else pure $ S.namedE "ERROR" S.SMat -- this shouldn't show up in stan code at all
  pure $ ProjModelData projData nNullVecsE' nPredictorsE' dmE pmE

-- given K null vectors, S states, and D predictors
-- alpha, theta, sigma
-- alpha is a K row-vector or S x K matrix
data Alpha = SimpleAlpha (S.Parameter S.ERVec) | HierarchicalAlpha (S.Parameter S.EMat)
-- theta is a D x K matrix (or Nothing)
newtype Theta = Theta (Maybe (S.Parameter S.EMat))
-- sigma is a K row-vector
newtype Sigma = Sigma {unSigma :: S.Parameter S.ERVec }

newtype Nu = Nu { unNu :: S.Parameter S.ERVec }

data ProjModelParameters where
  NormalProjModelParameters :: Alpha -> Theta -> Sigma -> ProjModelParameters
  CauchyProjModelParameters :: Alpha -> Theta -> Sigma -> ProjModelParameters
  StudentTProjModelParameters :: Alpha -> Theta -> Sigma -> Nu -> ProjModelParameters

paramTheta :: ProjModelParameters -> Theta
paramTheta (NormalProjModelParameters _ t _) = t
paramTheta (CauchyProjModelParameters _ t _) = t
paramTheta (StudentTProjModelParameters _ t _ _) = t

projModelParameters :: ModelConfig k md -> ProjModelData outerK md -> S.StanModelBuilderEff (ProjData outerK md) () ProjModelParameters
projModelParameters mc pmd = do
  let stdNormalDWA :: (S.TypeOneOf t [S.EReal, S.ECVec, S.ERVec], S.GenSType t) => S.DensityWithArgs t
      stdNormalDWA = S.DensityWithArgs S.std_normal TNil --(S.realE 0 :> S.realE 1 :> TNil)
      numPredictors = SBB.rowLength mc.designMatrixRow
  -- for now all the thetas are iid std normals

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
  let nStatesE = S.groupSizeE stateG
      hierAlphaNDS = S.NamedDeclSpec "alpha" $ S.matrixSpec nStatesE pmd.nNullVecsE
      fstI x k = S.sliceE S.s0 k x
      loopNVs = S.for "k" (S.SpecificNumbered (S.intE 1) pmd.nNullVecsE)
      diagPostMult m cv = S.diag_post_multiply m cv
      rowsOf nRowsE rv = S.diag_post_multiply (S.rep_matrix (S.realE 1) nRowsE (S.size rv))  $ S.transposeE rv
--      colsOf nColsE cv = diagPostMult (S.functionE S.rep_matrix (S.realE 1 :> S.functionE S.size (cv :> TNil) :> nColsE) cv :> TNil)
      hierAlphaPs = do
        muAlphaP <- S.simpleParameterWA
                    (S.NamedDeclSpec "muAlpha" $ S.rowVectorSpec pmd.nNullVecsE)
                    stdNormalDWA
        sigmaAlphaP <-  S.simpleParameterWA
                        (S.NamedDeclSpec "sigmaAlpha" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) $ S.rowVectorSpec pmd.nNullVecsE)
                        stdNormalDWA
        pure (muAlphaP :> sigmaAlphaP :> TNil)
  alpha <- case mc.alphaModel of
    AlphaSimple -> do
      fmap SimpleAlpha
        $ S.simpleParameterWA
           (S.NamedDeclSpec "alpha" $ S.rowVectorSpec pmd.nNullVecsE)
           stdNormalDWA
    AlphaHierCentered -> do
      alphaPs <- hierAlphaPs
      fmap HierarchicalAlpha
        $ S.addBuildParameter
        $ S.UntransformedP hierAlphaNDS [] alphaPs
        $ \(muAlphaE :> sigmaAlphaE :> TNil) m
          -> S.addStmt
             $ loopNVs
             $ \k -> S.sample (m `fstI` k) S.normalS (muAlphaE `fstI` k :> sigmaAlphaE `fstI` k :> TNil)
    AlphaHierNonCentered -> do
      alphaPs <- hierAlphaPs
      fmap HierarchicalAlpha
        $ S.withIIDRawMatrix hierAlphaNDS S.TransformedParametersBlock Nothing stdNormalDWA alphaPs
        $ \(muAlphaE :> sigmaAlphaE :> TNil) rawM -> rowsOf nStatesE muAlphaE `S.plusE` diagPostMult rawM (S.transposeE sigmaAlphaE)
  case mc.distribution of
    NormalDist -> pure $ NormalProjModelParameters alpha theta sigma
    CauchyDist -> pure $ CauchyProjModelParameters alpha theta sigma
    StudentTDist -> do
      let kVectorOf x = S.rep_row_vector (S.realE x) pmd.nNullVecsE
      nu <-  fmap Nu
             $ S.simpleParameterWA
             (S.NamedDeclSpec "nu" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) $ S.rowVectorSpec pmd.nNullVecsE)
             (S.DensityWithArgs S.gamma (kVectorOf 2 :> kVectorOf 0.1 :> TNil))
      pure $ StudentTProjModelParameters alpha theta sigma nu

data RunConfig = RunConfig { rcIncludePPCheck :: Bool, rcIncludeLL :: Bool }

-- not returning anything for now
projModel :: (Typeable outerK, Typeable md)
          => RunConfig -> ModelConfig k md -> ProjDataRTT outerK md -> S.StanModelBuilderEff (ProjData outerK md) () ()
projModel rc mc projData = do
  mData <- projModelData mc projData
  mParams <- projModelParameters mc mData
  let pExpr = S.parameterExpr
  let betaNDS = S.NamedDeclSpec "beta" $ S.matrixSpec mData.nPredictorsE mData.nNullVecsE
      nRowsE = S.dataSetSizeE mData.projDataTag
      fstI x k = S.sliceE S.s0 k x
      sndI x k = S.sliceE S.s1 k x
      loopNVs = S.for "k" (S.SpecificNumbered (S.intE 1) mData.nNullVecsE)
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
        $ \k -> let colk :: S.UExpr t -> S.UExpr (S.Sliced S.N1 t)
                    colk = flip sndI k --S.sliceE S.s1 k x
                in
                  S.grouped [ (sds `fstI` k) S.|=| S.sd (colk mData.projectionsE)
                            , colk stdNVPs S.|=| (colk mData.projectionsE |/| (sds `fstI` k))]
      let inverse :: (t ~ S.BinaryResultT S.BMultiply S.EReal t) => S.IntE -> S.UExpr t -> S.UExpr t --S.UExpr (TEO.BinaryResultT TEO.BMultiply S.EReal t)
          inverse k psCol = sds `fstI` k `S.timesE` psCol
      pure (stdNVPs, inverse)
    False -> pure (mData.projectionsE, const id)

  -- model
  let reIndexByState = S.indexE S.s0 (S.dataByGroupIndexE mData.projDataTag stateG)
      muE :: Alpha -> Theta -> S.IntE -> S.VectorE
      muE a t k =  case a of
       SimpleAlpha alphaP -> case t of
         Theta Nothing -> S.rep_vector (pExpr alphaP `fstI` k) nRowsE
         Theta (Just thetaP) -> pExpr alphaP `fstI` k `S.plusE` (predM `S.timesE` (pExpr thetaP `sndI` k))
       HierarchicalAlpha alphaP -> case t of
         Theta Nothing -> reIndexByState (pExpr alphaP `sndI` k)
         Theta (Just thetaP) -> reIndexByState (pExpr alphaP `sndI` k) `S.plusE` (predM `S.timesE` (pExpr thetaP `sndI` k))
      sigmaE :: Sigma -> S.IntE -> S.VectorE
      sigmaE s k = S.rep_vector (pExpr (unSigma s) `fstI` k) nRowsE

  let ppF :: Int
          -> ((S.IntE -> S.ExprList xs) -> S.IntE -> S.UExpr S.EReal)
          -> (S.IntE -> S.CodeWriter (S.IntE -> S.ExprList xs))
          -> S.StanModelBuilderEff (ProjData outerK md) () (S.ArrayE S.EReal)
      ppF k rngF rngPSCW =  S.generatePosteriorPrediction'
                            mData.projDataTag
                            (S.NamedDeclSpec ("predProj_" <> show k) $ S.array1Spec nRowsE S.realSpec)
                            rngF
                            (rngPSCW (S.intE k))
                            --               (pure $ \nE -> muE kE `fstI` nE :> unSigma mParams.pSigma `fstI` kE :> TNil)
                            (\_ p -> inverseF (S.intE k) p)
  let (sampleStmtF, ppStmtF) = case mParams of
        NormalProjModelParameters a t s ->
          let ssF e k = S.sample e S.normal (muE a t k :> sigmaE s k :> TNil)
              rF f nE = S.functionE S.normal_rngF (f nE)
              rpF k = pure $ \nE -> muE a t k `fstI` nE :> sigmaE s k `fstI` nE :> TNil
          in (ssF, \n -> ppF n rF rpF)
        CauchyProjModelParameters a t s ->
          let ssF e k = S.sample e S.cauchy (muE a t k :> sigmaE s k :> TNil)
              rF f nE = S.functionE S.cauchy_rngF (f nE)
              rpF k = pure $ \nE -> muE a t k `fstI` nE :> sigmaE s k `fstI` nE :> TNil
          in (ssF, \n -> ppF n rF rpF)
        StudentTProjModelParameters a t s n ->
          let nu :: Nu -> S.IntE -> S.VectorE
              nu n' k = S.rep_vector (pExpr (unNu n') `fstI` k) nRowsE
              ssF e k = S.sample e S.student_t (nu n k :> muE a t k :> sigmaE s k :> TNil)
              rF f nE = S.functionE S.student_t_rngF (f nE)
              rpF k  = pure $ \nE -> nu n k `fstI` nE :> muE a t k `fstI` nE :> sigmaE s k `fstI` nE :>  TNil
          in (ssF, \n' -> ppF n' rF rpF)

  S.inBlock S.SBModel $ S.addFromCodeWriter $ do
    let loopBody k = S.cwStmt_ $ S.addStmt $ sampleStmtF (nvps `sndI` k) k
    S.addStmt $ loopNVs loopBody
  -- generated quantities
  when rc.rcIncludePPCheck $ forM_ [1..modelNumNullVecs mc] ppStmtF
  pure ()

runProjModel :: forall (ks :: [(Symbol, Type)]) md r .
                (K.KnitEffects r
                , BRCC.CacheEffects r
                , ks F.⊆ DDP.ACSa5ByPUMAR
                , Typeable md
                , Flat.Flat (md [(Double, Double)])
                )
             => Bool
             -> RunConfig
             -> ModelConfig (F.Record ks) md
             -> DMS.MarginalStructure (Sum Double) (F.Record ks)
             -> FL.Fold (F.Record DDP.ACSa5ByPUMAR) (md Double)
             -> K.Sem r (K.ActionWithCacheTime r (ModelResult Text md)) -- no returned result for now
runProjModel clearCaches rc mc ms datFld = do
  let cacheDirE = (if clearCaches then Left else Right) "model/demographic/nullVecProjModel1_A5/"
      dataName = "projectionData_" <> dataText mc
  stanDir <- K.liftKnit MST.stanDir >>= K.knitMaybe "runModel: empty stanDir!" . BRCC.insureFinalSlash
  let runnerInputNames = S.RunnerInputNames
                         (stanDir <> "demographic/nullVecProj_M1_A5")
                         (modelText mc)
                         (Just $ S.GQNames "pp" dataName) -- posterior prediction vars to wrap
                         dataName
      (srcWindow, cachedSrc) = ACS.acs1Yr2012_21
  acsByPUMA_C <- DDP.cachedACSa5ByPUMA srcWindow cachedSrc 2021
--  acsByPUMA <- K.ignoreCacheTime acsByPUMA_C
  let outerKey = F.rcast @[GT.StateAbbreviation, GT.PUMA]
      catKey = F.rcast @ks
      count = realToFrac . view DT.popCount
      projData acsByPUMA =
        ProjData (modelNumNullVecs mc) (SBB.rowLength mc.designMatrixRow)
        (FL.fold (nullVecProjectionsModelDataFld ms mc.projVecs outerKey catKey count datFld) acsByPUMA)
  let modelData_C = fmap projData acsByPUMA_C
      meanSDFld :: FL.Fold Double (Double, Double) = (,) <$> FL.mean <*> FL.std
      meanSDFlds :: Int -> FL.Fold [Double] [(Double, Double)]
      meanSDFlds m = traverse (\n -> FL.premap (List.!! n) meanSDFld) [0..(m - 1)]
--        foldl' (\flds fld -> g <$> flds <*> fld) (pure []) $ replicate n meanSDFld
--        where g ls l = ls ++ [l]
  modelData <- K.ignoreCacheTime modelData_C
  let meanSDs = FL.fold (FL.premap (\(_, _, v) -> VS.toList v) $ meanSDFlds (modelNumNullVecs mc)) $ pdRows modelData
  K.logLE K.Info $ "meanSDs=" <> show meanSDs
  states <- FL.fold (FL.premap (view GT.stateAbbreviation) FL.set) <$> K.ignoreCacheTime acsByPUMA_C
  (dw, code) <-  S.dataWranglerAndCode modelData_C (pure ())
                 (stateGroupBuilder (view GT.stateAbbreviation)  (Set.toList states))
                 (const $ pure ())
                 (\projData' _ -> projModel rc mc projData')

  let nNullVecs = modelNumNullVecs mc
      unwraps = (\n -> S.UnwrapExpr ("matrix(ncol="
                                       <> show nNullVecs
                                       <> ", byrow=TRUE, unlist(jsonData $ nvp_ProjectionData))[,"
                                       <> show n <> "]") ("obsNVP_" <> show n))
                <$> [1..nNullVecs]
  res_C <- S.runModel' @BRCC.SerializerC @BRCC.CacheData
           cacheDirE
           (Right runnerInputNames)
           Nothing
           dw
           code
           (projModelResultAction mc) --S.DoNothing -- (stateModelResultAction mcWithId dmr)
           (S.ShinyStan unwraps) --(S.Both [S.UnwrapNamed "successes" "yObserved"])
           modelData_C
           (pure ())
  K.logLE K.Info "projModel run complete."
  pure res_C

--NB: parsed summary data has stan indexing, i.e., Arrays start at 1.
projModelResultAction :: forall outerK md k r .
                         (K.KnitEffects r
                         , Typeable md
                         , Typeable outerK
                         )
                      => ModelConfig k md
                      -> S.ResultAction (ProjData outerK md) () S.DataSetGroupIntMaps S.DataSetGroupIntMaps r () (ModelResult Text md)
projModelResultAction mc = S.UseSummary f where
  f summary _ modelDataAndIndexes_C _ = do
    (modelData, resultIndexesE) <- K.ignoreCacheTime modelDataAndIndexes_C
    -- compute means of predictors because model was zero-centered in them
    let nPredictors = SBB.rowLength mc.designMatrixRow
        mdMeansFld = FL.premap (\(_, md, _) -> mc.mdFuncs.mdfToList md)
                    $ traverse (\n -> FL.premap (List.!! n) FL.mean) [0..(nPredictors - 1)]
        nvpSDFld = FL.premap (\(_, _, v) -> VS.toList v)
                   $ traverse (\n -> FL.premap (List.!! n) FL.std) [0..((modelNumNullVecs mc) - 1)]
        (mdMeansL, nvpSDsL) = FL.fold ((,) <$> mdMeansFld <*> nvpSDFld) $ pdRows modelData
        rescaleAlphaBeta xs = if mc.standardizeNVs then zipWith (*) xs nvpSDsL else xs
    stateIM <- K.knitEither
      $ resultIndexesE >>= S.getGroupIndex (S.RowTypeTag @(ProjDataRow outerK md) "ProjectionData") stateG
    let allStates = IM.elems stateIM
        getVector n = K.knitEither $ S.getVector . fmap CS.mean <$> S.parse1D n (CS.paramStats summary)
        getMatrix n = K.knitEither $ fmap CS.mean <$> S.parse2D n (CS.paramStats summary)
    geoMap <- case mc.alphaModel of
      AlphaSimple -> do
        alphaV <- getVector "alpha" -- states by nNullvecs
        pure $ M.fromList $ fmap (, rescaleAlphaBeta $ V.toList alphaV) allStates
      _ -> do
        alphaVs <- getMatrix "alpha"
        let mRowToList cols row = fmap (\c -> S.getIndexed alphaVs (row, c)) $ [1..cols]
        pure $ M.fromList $ fmap (\(row, sa) -> (sa, rescaleAlphaBeta $ mRowToList (modelNumNullVecs mc) row)) $ IM.toList stateIM
    betaSIL <- case nPredictors of
      0 -> pure $ replicate (modelNumNullVecs mc) []
      p -> do
        betaVs <- getMatrix "beta" -- ps by nNullVecs
        let mColToList rows col = fmap (\r -> S.getIndexed betaVs (r, col)) [1..rows]
        pure $ transp $ fmap (\m -> zip (rescaleAlphaBeta $ mColToList p m) mdMeansL) [1..modelNumNullVecs mc]
    betaSI <- K.knitEither $ mc.mdFuncs.mdfFromList betaSIL
    pure $ ModelResult geoMap betaSI

transp :: [[a]] -> [[a]]
transp = go [] where
  go x [] = fmap reverse x
  go [] (r : rs) = go (fmap (: []) r) rs
  go x (r :rs) = go (List.zipWith (:) r x) rs


data ASERModelP a = ASERModelP { mASER_PWLogDensity :: a, mASER_FracOver45 :: a, mASER_FracGrad :: a, mASER_FracOfColor :: a , mASER_FracWNG :: a  }
  deriving stock (Show, Generic)
  deriving anyclass Flat.Flat

a4serModelDatFld :: (F.ElemOf rs DT.PWPopPerSqMile
                    , F.ElemOf rs DT.Age4C
                    , F.ElemOf rs DT.Education4C
                    , F.ElemOf rs DT.Race5C
                    , F.ElemOf rs DT.PopCount
                   )
                 => FL.Fold (F.Record rs) (ASERModelP Double)
a4serModelDatFld = ASERModelP <$> dFld <*> aFld <*> gFld <*> rFld <*> wngFld
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

a5serModelDatFld :: (F.ElemOf rs DT.PWPopPerSqMile
                    , F.ElemOf rs DT.Age5FC
                    , F.ElemOf rs DT.Education4C
                    , F.ElemOf rs DT.Race5C
                    , F.ElemOf rs DT.PopCount
                    )
                 => FL.Fold (F.Record rs) (ASERModelP Double)
a5serModelDatFld = ASERModelP <$> fmap safeLog dFld <*> aFld <*> gFld <*> rFld <*> wngFld
  where
    nPeople = realToFrac . view DT.popCount
    safeLog x = if x > 1 then Numeric.log x else 0
    dens = view DT.pWPopPerSqMile --let x = view DT.pWPopPerSqMile r in
    wgtFld = FL.premap nPeople FL.sum
    wgtdFld f = safeDiv <$> FL.premap (\r -> nPeople r * f r) FL.sum <*> wgtFld
    dFld = wgtdFld dens
    over45 = (`elem` [DT.A5F_45To64, DT.A5F_65AndOver]) . view DT.age5FC
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
    fracFld f = safeDiv <$> FL.prefilter f wgtFld <*> wgtFld
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
