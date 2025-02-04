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
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE UnicodeSyntax #-}
{-# LANGUAGE StandaloneDeriving #-}

module BlueRipple.Model.Demographic.BLCorrModel
  (
    module BlueRipple.Model.Demographic.BLCorrModel
  )
where

import qualified BlueRipple.Data.CachingCore as BRCC
import qualified BlueRipple.Model.Demographic.DataPrep as DDP
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

import qualified Frames as F
import qualified Frames.Melt as F
import qualified Frames.Serialize as FS
import qualified Frames.Transform as FT
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vinyl as V
import qualified Data.Vinyl.TypeLevel as V

import Control.Lens (view)

import qualified Stan as S
import qualified Stan.BuildingBlocks as SBB (rowLength)
import Stan (TypedList(..))
import Stan.Operators

import qualified Flat

data PopAndDensity = PopAndDensity { pop :: Int, pwDensity :: Double}
instance Semigroup PopAndDensity where
  (<>) (PopAndDensity pa pwda) (PopAndDensity pb pwdb) = PopAndDensity p pwd where
    p = pa + pb
    pwd = if p > 0 then (realToFrac pa * pwda + realToFrac pb * pwdb) / realToFrac p else 0

instance Monoid PopAndDensity where
  mempty = PopAndDensity 0 0
  mappend = (<>)

marginalFld :: (Ord ck, Ord mk, Monoid d)
            => (F.Record qs -> ck)
            -> (F.Record qs -> mk)
            -> (F.Record qs -> d)
            -> FL.Fold (F.Record qs) [(ck, Map mk d)]
marginalFld catKeyF mKeyF datF =
  MR.mapReduceFold
  MR.noUnpack
  (MR.assign catKeyF (\r -> (mKeyF r, datF r)))
  (MR.foldAndLabel (FL.foldByKeyMap FL.mconcat) (,))

data DataRow rs = DataRow (F.Record rs) [Int]
deriving instance (Show (F.Record rs)) => Show (DataRow rs)

--type instance S.DataSource S.ModelDataT = [DataRows] --F.Frame FB_Result
--type instance S.DataSource S.GQDataT = () --F.Frame FB_Matchup

makeRowFromPD :: F.Record cs -> [PopAndDensity] -> DataRow (cs V.++ '[DT.PWPopPerSqMile])
makeRowFromPD catR dats = DataRow (catR F.<+> FT.recordSingleton @DT.PWPopPerSqMile pwd) (fmap pop dats) where
  pwd = pwDensity $ mconcat dats

dataRowsFld :: forall rs ck k qs d .
               (Ord (F.Record rs)
               , Ord ck
               , Ord k
               , Monoid d
               , BRK.FiniteSet k
               )
            => (F.Record qs -> ck)
            -> (F.Record qs -> k)
            -> (F.Record qs -> d)
            -> (ck -> [d] -> DataRow rs)
            -> FL.Fold (F.Record qs) [DataRow rs]
dataRowsFld catKeyF mKeyF datF mkRow = fmap (\(qs, m) -> mkRow qs (M.elems $ M.union m zeroMap)) <$> marginalFld catKeyF mKeyF datF
  where
    zeroMap = M.fromList $ fmap (,mempty) (S.toList $ BRK.elements @k)

dataRowRec :: DataRow rs -> F.Record rs
dataRowRec (DataRow r _) = r

dataRowCounts :: DataRow rs -> [Int]
dataRowCounts (DataRow _ ms) = ms

instance (V.RMap rs, FS.RecFlat rs) => Flat.Flat (DataRow rs) where
  size (DataRow r ms) = Flat.size (FS.toS r, ms)
  encode (DataRow r ms) = Flat.encode (FS.toS r, ms)
  decode = fmap (\(sr, ms) -> DataRow (FS.fromS sr) ms) Flat.decode

type DataRows rs = [DataRow rs]
type DataRTT rs = S.RowTypeTag (DataRow rs)

modelIDT :: S.InputDataType S.ModelDataT (DataRows rs)
modelIDT = S.ModelData

--deriving anyclass instance (Flat.Flat (ProjDataRow rs)) => Flat.Flat (ProjData rs)
{-
--Figure out results section once we see if this model is well-behaved

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
-}

data CovarianceStructure = DiagonalCovariance | LKJCovariance Int deriving stock (Show, Eq, Generic)

covarianceText :: CovarianceStructure -> Text
covarianceText DiagonalCovariance = "diagCov"
covarianceText (LKJCovariance n) = "lkj" <> show n

data BetaModel = BetaSimple | BetaHierCentered CovarianceStructure | BetaHierNonCentered CovarianceStructure deriving stock (Show, Eq, Generic)

betaModelText :: BetaModel -> Text
betaModelText BetaSimple = "BS"
betaModelText (BetaHierCentered cs) = "BHC_" <> covarianceText cs
betaModelText (BetaHierNonCentered cs) = "BHNC_" <> covarianceText cs

data ModelConfig alphaK (pd :: Type -> Type) where
  ModelConfig :: ({-Traversable pd-})
              => { nCounts :: Int
                 , alphaDMR :: S.DesignMatrixRow alphaK
--                 , predDMR :: DM.DesignMatrixRow (pd Double)
                 , betaModel :: BetaModel
                 , betaLastAsZero :: Bool
                 , dirichletPrior :: Bool
                 } -> ModelConfig alphaK pd

modelText :: ModelConfig alphaK md -> Text
modelText mc = mc.alphaDMR.dmName
--               <> "_" <> mc.predDMR.dmName
               <> "_" <> betaModelText mc.betaModel
               <> if mc.betaLastAsZero then "_l0" else ""
               <> if mc.dirichletPrior then "_dir" else ""

dataText :: ModelConfig alphaK md -> Text
dataText mc = mc.alphaDMR.dmName
--              <> "_" <> mc.predDMR.dmName

stateG :: S.GroupTypeTag Text
stateG = S.GroupTypeTag "State"

stateGroupBuilder :: forall f rs. (Foldable f, Typeable rs)
                  => (F.Record rs -> Text) -> f Text -> S.StanDataBuilderEff S.ModelDataT (DataRows rs) (DataRTT rs)
stateGroupBuilder saF states = do
  dataSetTag <- S.addData "CountData" (modelIDT @rs) (S.ToFoldable id)
  S.addGroupIndexForData (modelIDT @rs) stateG dataSetTag $ S.makeIndexFromFoldable show (saF . dataRowRec) states
  S.addGroupIntMapForData (modelIDT @rs) stateG dataSetTag $ S.dataToIntMapFromFoldable (saF . dataRowRec) states
  pure dataSetTag

data ModelData rs =
  ModelData
  {
    dataTag :: S.RowTypeTag (DataRow rs)
  , nCatsE :: S.IntE
  , countsE :: S.ArrayE (S.EArray1 S.EInt)
  , nCovariatesE :: S.IntE
  , covariatesE :: S.MatrixE
--  , nPredictorsE :: TE.IntE
--  , predictorsE :: TE.MatrixE
  }

--TODO: add predictors to alphas to make one matrix of covariates
modelData :: forall pd alphaK rs . (Typeable rs)
          => ModelConfig alphaK pd
          -> S.RowTypeTag (DataRow rs)
          -> (F.Record rs -> alphaK)
          -> (F.Record rs -> pd Double)
          -> S.StanModelBuilderEff (DataRows rs) () (ModelData rs)
modelData mc dataSetTag catKey _predF = do
--  dat <- S.dataSetTag @(DataRow rs) SC.ModelData "CountData"
  (countsE', nCatsE') <- S.addArrayOfIntArrays (modelIDT @rs) dataSetTag "MCounts" Nothing mc.nCounts dataRowCounts (Just 0) Nothing
  let (_, nCovariatesE') = S.designMatrixColDimBinding mc.alphaDMR Nothing
  covariatesDME <- if SBB.rowLength mc.alphaDMR > 0
                   then S.addDesignMatrix (modelIDT @rs) dataSetTag (contramap (catKey . dataRowRec) mc.alphaDMR) Nothing
                   else pure $ S.namedE "ERROR" S.SMat -- this shouldn't show up in stan code at all
{-  let (_, nPredictorsE') = DM.designMatrixColDimBinding mc.predDMR Nothing
  dmE <- if DM.rowLength mc.predDMR > 0
         then DM.addDesignMatrix dat (contramap (predF . projRowRec) mc.predDMR) Nothing
         else pure $ TE.namedE "ERROR" TE.SMat -- this shouldn't show up in stan code at all
-}
  pure $ ModelData dataSetTag nCatsE' countsE' nCovariatesE' covariatesDME --nPredictorsE' dmE

-- S states
-- K categories to count
-- C categories to use for prediction
-- D predictors. 0 for now.
-- M is number of one-hot encoded alphas
-- M x K matrix or S array of M x K matrix. M=Number of cols is < 1 + C + D since we binary or one-hot encode all the categories
data Beta = SimpleBeta (S.Parameter S.EMat) | HierarchicalBeta (S.Parameter (S.EArray1 S.EMat))

data ModelParameters where
  ModelParameters :: Beta -> ModelParameters

modelBeta :: ModelConfig alphaK pd -> ModelData rs -> S.StanModelBuilderEff (DataRows rs) () Beta
modelBeta mc pmd = do
  let nStatesE = S.groupSizeE stateG
--      toVector x = S.to_vector (x :> TNil)
      betaColsE = if mc.betaLastAsZero then pmd.nCatsE |-| S.intE 1 else pmd.nCatsE
      betaName = if mc.betaLastAsZero then "betaR" else "beta"
      betaShape = S.matrixSpec pmd.nCovariatesE betaColsE
      hierBetaSpec =  S.array1Spec nStatesE betaShape
      hierBetaPs :: S.StanModelBuilderEff (DataRows rs) () (S.Parameters [S.EMat, S.EMat])
      hierBetaPs = do
        muBetaP <- S.addBuildParameter
                    $ S.UntransformedP
                    (S.NamedDeclSpec ("mu" <> betaName) betaShape)
                    [] TNil
                    (\_ p -> S.addStmt $ S.sample (S.to_vector p) S.std_normal TNil)

        tauBetaP <- S.addBuildParameter
                     $ S.UntransformedP
                     (S.NamedDeclSpec ("tau" <> betaName) $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) betaShape)
                     [] TNil
                     (\_ p -> S.addStmt $ S.sample (S.to_vector p) S.std_normal TNil)
        pure (muBetaP :> tauBetaP :> TNil)
      betaHier cs cent = do
        hierPs <- hierBetaPs
        case hierPs of
          (muBetaP :> tauBetaP :> TNil) -> do
            muBetaAP <- S.addBuildParameter
                       $ S.TransformedP
                       (S.NamedDeclSpec ("mu" <> betaName <> "A") $ S.array1Spec nStatesE betaShape)
                       []
                       (muBetaP :> TNil)
                       S.TransformedParametersBlock
                       (\(muBetaE :> TNil) -> S.DeclRHS $ S.rep_array1 muBetaE nStatesE)
                       TNil
                       (\_ _ -> pure ())

{-
              SMB.addFromCodeWriter
                        $ TE.declareRHSNW (TE.NamedDeclSpec ("mu" <> betaName <> "A") $ TE.array1Spec nStatesE $ betaShape [])
                        $ TE.functionE SF.rep_array (DAG.parameterExpr muBetaP :> (nStatesE :> TNil))
-}
--            let tauBeta = DAG.parameterExpr tauBetaP
            case cs of
              DiagonalCovariance -> do
                fmap HierarchicalBeta
                  $ S.matrixMultiNormalParameter' S.Diagonal cent muBetaAP tauBetaP
                  (S.NamedDeclSpec betaName $ hierBetaSpec)
              LKJCovariance lkjPriorP -> do
                lkjCorrBetaP <- S.simpleParameter
                                (S.NamedDeclSpec ("lkj" <> betaName)
                                 $ S.choleskyFactorCorrSpec (pmd.nCovariatesE |*| betaColsE))
                                (S.given (S.realE $ realToFrac lkjPriorP) :> TNil)
                                S.lkj_corr_cholesky
                fmap HierarchicalBeta
                  $ S.matrixMultiNormalParameter' (S.Cholesky lkjCorrBetaP) cent muBetaAP tauBetaP
                  (S.NamedDeclSpec betaName $ hierBetaSpec)
--          _ -> SMB.stanBuildError "BLCorrModel.modelBeta: Pattern match error in hierarchical beta paramters. Yikes."
  betaRawP <- case mc.betaModel of
    BetaSimple -> fmap SimpleBeta
                  $ S.addBuildParameter
                  $ S.UntransformedP (S.NamedDeclSpec betaName $ S.matrixSpec pmd.nCovariatesE betaColsE)
                  [] TNil (\_ muM -> S.addStmt $ S.sample (S.to_vector muM) S.std_normal TNil)
    BetaHierCentered cs -> betaHier cs S.Centered
    BetaHierNonCentered cs -> betaHier cs S.NonCentered
  case  mc.betaLastAsZero of
    False -> pure betaRawP
    True -> do
      let fullBetaShape = S.matrixSpec pmd.nCovariatesE pmd.nCatsE
          zeroCol = S.rep_vector (S.realE 0) pmd.nCovariatesE
          appendZeroCol m = S.append_col m zeroCol
      case betaRawP of
        SimpleBeta brP ->
          fmap SimpleBeta
          $ S.addBuildParameter
          $ S.TransformedP (S.NamedDeclSpec "beta" $ fullBetaShape) []
          (brP :> TNil) S.TransformedParametersBlock
          (\(br :> TNil) -> S.DeclRHS $ appendZeroCol br)
          TNil
          (\_ _ -> pure ())
        HierarchicalBeta brP ->
          fmap HierarchicalBeta
          $ S.addBuildParameter
          $ S.TransformedP (S.NamedDeclSpec "beta" $ S.array1Spec nStatesE $ fullBetaShape) []
          (brP :> TNil) S.TransformedParametersBlock
          (\(br :> TNil) -> S.DeclCodeF
            $ \b -> S.addStmt $ S.loopSized nStatesE "s"
                    $ \ns -> b !! ns |=| appendZeroCol (br !! ns))
          TNil
          (\_ _ -> pure ())

data RunConfig = RunConfig { rcIncludePPCheck :: Maybe Int, rcIncludeLL :: Bool, statesM :: Maybe (Text, [Text]) }

projModel :: Typeable rs
          => RunConfig
          -> (F.Record rs -> alphaK)
          -> (F.Record rs -> pd Double)
          -> ModelConfig alphaK pd
          -> DataRTT rs
          -> S.StanModelBuilderEff  (DataRows rs) () ()
projModel rc alphaKeyF predF mc dataRtt = do
  mData <- modelData mc dataRtt alphaKeyF predF
  let nRowsE = S.dataSetSizeE mData.dataTag
  -- transformed data
  totalCountE <- S.inBlock S.SBTransformedDataGQ $ S.addFromCodeWriter $ do
    tc <- S.declareNW
      (S.NamedDeclSpec "TCount" $ S.array1Spec nRowsE $ S.addVMs (S.Modifiers [S.lowerM $ S.intE 0]) S.intSpec)
    S.addStmt $ S.loopSized nRowsE "n" $ \n -> (tc !! n) |=| S.sumInt (mData.countsE !! n)
    pure tc
  betaP <- modelBeta mc mData
  let reIndexByState = S.indexE S.s0 (S.dataByGroupIndexE mData.dataTag stateG)
      betaByRow :: S.IntE -> S.MatrixE
      betaByRow ie = case betaP of
          SimpleBeta m ->  S.parameterExpr m
          HierarchicalBeta betaByState -> reIndexByState (S.parameterExpr betaByState) !! ie
      mnArgByRowE nE = S.transposeE $ (mData.covariatesE !! nE) |*| betaByRow nE


  case mc.dirichletPrior of
    False -> do
      S.inBlock S.SBModel
        $ S.addFromCodeWriter
        $ S.addStmt
        $ S.loopSized nRowsE "n"
        $ \nE -> S.target $ S.densityE S.multinomial_logit_lpmf (mData.countsE !! nE) (mnArgByRowE nE :> TNil)

      case rc.rcIncludePPCheck of
        Just nChoices -> do
          let ppCheck n = S.inBlock S.SBGeneratedQuantities
                          $ S.generatePosteriorPrediction'
                          mData.dataTag
                          (S.NamedDeclSpec "ppCounts"
                           $ S.array1Spec nRowsE S.intSpec
                          )
                          (\f nE -> (S.functionE S.multinomial_logit_rngF (mnArgByRowE nE :> f nE)) !! S.intE n)
                          (pure $ \nE -> totalCountE !! nE :> TNil)
                          (const id)
          mapM_ ppCheck [1..nChoices]
        Nothing -> pure ()
      when rc.rcIncludeLL
        $ S.generateLogLikelihood
        mData.dataTag
        S.multinomialLogitDist
        (pure $ \nE -> mnArgByRowE nE :> TNil)
        (pure $ \nE -> mData.countsE !! nE)

    True -> do
      (dirichlet_multinomial, dirichlet_multinomial_lpmf, dirichlet_multinomial_rngF, _) <- S.dirichletMultinomial @_ @S.ECVec
--      let softmax x = TE.functionE SF.softmax (x :> TNil)
      dPrecE <- fmap S.parameterExpr
                $ S.addBuildParameter
                $ S.UntransformedP
                (S.NamedDeclSpec "dPrec" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)  [] TNil
                (\_ dp -> S.addStmt $ S.sample dp S.normal (S.realE 5 :> S.realE 25 :> TNil))
      let dmArgByRowE nE = dPrecE |*| S.softmax (mnArgByRowE nE)
      S.inBlock S.SBModel
        $ S.addFromCodeWriter
        $ S.addStmt
        $ S.loopSized nRowsE "n"
        $ \nE -> S.cwStmt_
                 $ (do
                       S.addStmt $ S.sample (mData.countsE !! nE) dirichlet_multinomial (dmArgByRowE nE :> TNil)
                   )
      case  rc.rcIncludePPCheck of
        Just nChoices -> do
          let ppCheck n = S.inBlock S.SBGeneratedQuantities
                          $ S.generatePosteriorPrediction'
                          mData.dataTag
                          (S.NamedDeclSpec ("ppCounts_" <> show n)
                           $ S.array1Spec nRowsE S.intSpec
                          )
                          (\f nE -> (S.functionE dirichlet_multinomial_rngF (dmArgByRowE nE :> f nE)) !! S.intE n)
                          (pure $ \nE -> totalCountE !! nE :> TNil)
                          (const id)
          mapM_ ppCheck [1..nChoices]
        Nothing -> pure ()
      when rc.rcIncludeLL
        $ S.generateLogLikelihood'
        $ S.addToLLSet mData.dataTag
        (S.LLDetails
         (S.densityE dirichlet_multinomial_lpmf)
         (pure $ \nE -> dmArgByRowE nE :> TNil)
         (pure $ \nE -> mData.countsE !! nE)
        )
        S.emptyLLSet


-- rs is stateAbbr ++ kP ++ PWPopPerSqMile
runProjModel :: forall kM pd kPs rs r .
                (K.KnitEffects r
                , BRCC.CacheEffects r
--                , kP F.⊆ DDP.ACSByPUMAR
                , Typeable pd
                , Ord kM
                , BRK.FiniteSet kM
                , rs ~ kPs V.++ '[DT.PWPopPerSqMile]
                , V.RMap rs
                , FS.RecFlat rs
                , Ord (F.Record rs)
                , Show (F.Record rs)
                , Typeable rs
                , F.ElemOf rs GT.StateAbbreviation
                , Ord (F.Record kPs)
                , kPs F.⊆ DDP.ACSa5ByPUMAR
                )
             => Bool
             -> RunConfig
             -> ModelConfig (F.Record rs) pd
             -> (F.Record DDP.ACSa5ByPUMAR -> kM)
             -> (F.Record DDP.ACSa5ByPUMAR -> F.Record kPs)
             -> (F.Record rs -> pd Double)
             -> K.Sem r (K.ActionWithCacheTime r ())
runProjModel clearCaches rc mc margKeyF _predKeyF predF = do
  let cacheRoot = "model/demographic/nullVecProjModel/"
      cacheDirE = (if clearCaches then Left else Right) cacheRoot
      dataName = "blCorrData_" <> dataText mc <> maybe "" fst rc.statesM
      countF r = PopAndDensity (view DT.popCount r) (view DT.pWPopPerSqMile r)
  stanDir <- K.liftKnit MST.stanDir >>= K.knitMaybe "runModel: empty stanDir!" . BRCC.insureFinalSlash
  let runnerInputNames = S.RunnerInputNames
                         (stanDir <> "demographic/blCorrModel")
                         (modelText mc)
                         (Just $ S.GQNames "pp" dataName) -- posterior prediction vars to wrap
                         dataName
      statesFilter = maybe id (\(_, sts) -> F.filterFrame ((`elem` sts) . view GT.stateAbbreviation)) rc.statesM
      (srcWindow, cachedSrc) = ACS.acs1Yr2012_21
  acsByPUMA_C <- fmap statesFilter <$> DDP.cachedACSa5ByPUMA srcWindow cachedSrc 2021 -- most recent available
  let dataCacheKey = cacheRoot <> "/acsCounts_" <> mc.alphaDMR.dmName <> maybe "" fst rc.statesM
  when clearCaches $ BRCC.clearIfPresentD dataCacheKey
  acsCountedByPUMA_C <- BRCC.retrieveOrMakeD
                       dataCacheKey
                       acsByPUMA_C
                       $
                       \acsByPUMA -> do
                         let mkRow acsByPUMARow = makeRowFromPD (acsByPUMARow)
                             counted = FL.fold (dataRowsFld (F.rcast @kPs) margKeyF countF mkRow) acsByPUMA
--                         K.logLE K.Info $ "counted: " <> show counted
                         pure counted
  states <-  FL.fold (FL.premap (view GT.stateAbbreviation) FL.set) <$> K.ignoreCacheTime acsByPUMA_C
  (dw, code) <-  S.dataWranglerAndCode acsCountedByPUMA_C (pure ())
                 (stateGroupBuilder (view GT.stateAbbreviation)  (S.toList states))
                 (const $ pure ())
                 (\dataRTT _ -> projModel rc id predF mc dataRTT)

  let unwraps = case rc.rcIncludePPCheck of
        Just nChoices ->
          let f n = S.UnwrapExpr ("matrix(ncol="
                                   <> show nChoices
                                   <> ", byrow=TRUE, unlist(jsonData $ MCounts))[,"
                                   <> show n <> "]") ("yCounts_" <> show n)
          in fmap f [1..nChoices]
        Nothing -> []
--      unwraps = [SR.UnwrapNamed "MCounts" "yCounts"]
  res_C <- S.runModel' @BRCC.SerializerC @BRCC.CacheData
           cacheDirE
           (Right runnerInputNames)
           (Just $ S.StanMCParameters 4 4 (Just 1000) (Just 1000) Nothing Nothing (Just 1))
           dw
           code
           S.DoNothing
           (S.Both unwraps) --(SMR.Both [SR.UnwrapNamed "successes" "yObserved"])
           acsCountedByPUMA_C
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
--    eRP = DM.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)

designMatrixRow_E :: S.DesignMatrixRow (F.Record '[DT.Education4C])
designMatrixRow_E = S.DesignMatrixRow "E" [eRP]
  where
    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)


designMatrixRow_1_E :: S.DesignMatrixRow (F.Record '[DT.Education4C])
designMatrixRow_1_E = S.DesignMatrixRow "I_E" [cRP, eRP]
  where
    cRP = S.DesignMatrixRowPart "Ones" 1 (const $ VU.singleton 1) -- for pure (state-level) alpha
    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)

designMatrixRow_1_S_E :: S.DesignMatrixRow (F.Record '[DT.SexC, DT.Education4C])
designMatrixRow_1_S_E = S.DesignMatrixRow "I_S_E" [cRP, sRP, eRP]
  where
    cRP = S.DesignMatrixRowPart "Ones" 1 (const $ VU.singleton 1) -- for pure (state-level) alpha
    sRP = S.boundedEnumRowPart Nothing "Sex" (view DT.sexC)
    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)

designMatrixRow_1_S_E_R :: S.DesignMatrixRow (F.Record [DT.SexC, DT.Education4C, DT.Race5C])
designMatrixRow_1_S_E_R = S.DesignMatrixRow "I_S_E_R" [cRP, sRP, eRP, rRP]
  where
    cRP = S.DesignMatrixRowPart "Ones" 1 (const $ VU.singleton 1) -- for pure (state-level) alpha
    sRP = S.boundedEnumRowPart Nothing "Sex" (view DT.sexC)
    eRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Edu" (view DT.education4C)
    rRP = S.boundedEnumRowPart (Just DT.R5_WhiteNonHispanic) "Race" (view DT.race5C)
