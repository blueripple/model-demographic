{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE UnicodeSyntax #-}

module BlueRipple.Model.Demographic.StanModels
  (
    module BlueRipple.Model.Demographic.StanModels
  )
where

import qualified BlueRipple.Model.Demographic.DataPrep as DDP
import qualified BlueRipple.Data.Small.DataFrames as BRDF
import qualified BlueRipple.Data.Types.Demographic as DT
import qualified BlueRipple.Data.Types.Geographic as GT
import qualified BlueRipple.Data.ACS_PUMS as ACS
import qualified BlueRipple.Data.Keyed as BRK
import qualified BlueRipple.Data.CachingCore as BRKU

import qualified Stan as S
import Stan (TypedList(..))
import Stan.Operators
import qualified CmdStan as CS

--import qualified Frames.Streamly.InCore as FS
import qualified Frames.Serialize as FS
import qualified Flat
import Control.Lens (view, (^.))
import qualified Control.Foldl as FL
import qualified Data.Map as M
import qualified Data.Set as Set
import qualified Data.Vinyl as V
import qualified Data.Vinyl.TypeLevel as V
import qualified Data.Vector as Vec
import qualified Data.Vector.Unboxed as VU
import qualified Frames as F
import qualified Frames.Melt as F
import qualified Knit.Report as K
import qualified Numeric
import qualified Data.IntMap.Strict as IM

logLengthC :: (K.KnitEffects r, Foldable f) => K.ActionWithCacheTime r (f a) -> Text -> K.Sem r ()
logLengthC xC t = K.ignoreCacheTime xC >>= \x -> K.logLE K.Info $ t <> " has " <> show (FL.fold FL.length x) <> " rows."

runModel :: forall ks l r .
            (K.KnitEffects r, BRKU.CacheEffects r
            , Ord (F.Record ks)
            , Enum l, Bounded l, Ord l
            , Typeable (ks V.++ '[DT.PWPopPerSqMile])
            , F.ElemOf  (ks V.++ '[DT.PWPopPerSqMile]) DT.PWPopPerSqMile
            , ks F.⊆ ([BRDF.Year, GT.StateAbbreviation, GT.StateFIPS] V.++ (ks V.++ '[DT.PWPopPerSqMile]))
            , V.RMap ks
            , FS.RecFlat ks
            , Ord (F.Rec FS.SElField ks)
            , BRK.FiniteSet (F.Record ks)
            )
         => Bool
         -> ModelConfig ()
         -> (Text, F.Record DDP.ACSa6ByStateR -> l)
         -> (Text, F.Record DDP.ACSa6ByStateR -> F.Record ks, S.DesignMatrixRow (F.Record ks))
         -> K.Sem r (K.ActionWithCacheTime r (ModelResult Text ks))
runModel clearCaches mc (modeledT, modeledK) (fromT, cKey, dmr) = do
  let cacheDirE = let k = ("model/demographic/" <> modeledT <> "/") in if clearCaches then Left k else Right k
      dataName = "acs" <> modeledT <> "_" <> S.dmName dmr <> modelConfigSuffix mc
      runnerInputNames = S.RunnerInputNames
                         ("br-2022-Demographics/stan" <> modeledT)
                         ("normal" <> fromT <> "_" <> S.dmName dmr <> modelConfigSuffix mc)
                         (Just $ S.GQNames "pp" dataName)
                         dataName
      (srcWindow, cachedSrc) = ACS.acs1Yr2012_21
  acs_C <- DDP.cachedACSa6ByState srcWindow cachedSrc 2021 -- most recent available
--  K.ignoreCacheTime acs_C >>= BRK.logFrame
  logLengthC acs_C "acsByState"
  let acsMN_C = fmap (DDP.acsByStateMN cKey modeledK) acs_C
      mcWithId = "normal" <$ mc
--  K.ignoreCacheTime acsMN_C >>= print
  logLengthC acsMN_C ("acsByState Counted for " <> modeledT)
  states <- FL.fold (FL.premap (view GT.stateAbbreviation . fst) FL.set) <$> K.ignoreCacheTime acsMN_C
  (dw, code) <- S.dataWranglerAndCode acsMN_C (pure ())
                (groupBuilderState (Set.toList states))
                (const $ pure ()) -- no GQ group setup
                (\acsTag _ -> normalModel (contramap F.rcast dmr) mc acsTag)
  res_C <-S.runModel' @BRKU.SerializerC @BRKU.CacheData
          cacheDirE
          (Right runnerInputNames)
          Nothing
          dw
          code
          (stateModelResultAction mcWithId dmr)
          (S.Both [S.UnwrapNamed "successes" "yObserved"])
          acsMN_C
          (pure ())
  K.logLE K.Info "citizenModel run complete."
  pure res_C

data HierarchicalType = HCentered | HNonCentered deriving stock (Show, Eq, Ord)

data ModelConfig a = ModelConfig { modelID :: a
                                 , includeAlpha0 :: Bool
                                 , alphaType :: HierarchicalType
                                 , includeDensity :: Bool
                                 } deriving stock (Functor, Show)

modelConfigSuffix :: ModelConfig a -> Text
modelConfigSuffix (ModelConfig _ ia at id') = a0s <> ats <> ids
  where
    a0s = if ia then "_a0" else ""
    ats = if at == HCentered then "_ac" else "_anc"
    ids = if id' then "_d" else ""

modelName :: ModelConfig Text -> Text
modelName mc = modelID mc <> modelConfigSuffix mc

addTermMaybe :: Maybe a -> (a -> S.UExpr t -> S.UExpr t) -> S.UExpr t -> S.UExpr t
addTermMaybe mA combine e = case mA of
  Nothing -> e
  Just a -> combine a e

type Row rs a = (F.Record rs, a)
type ACSRowTag rs a = S.RowTypeTag (Row rs a)

modelIDT :: forall rs a . S.InputDataType S.ModelDataT [Row rs a]
modelIDT = S.ModelData

data ModelData rs = ModelData
  {
    acsDataTag :: ACSRowTag rs (VU.Vector Int)
  , nData :: S.IntE
  , nStates :: S.IntE
  , nPredictors :: S.IntE
  , trials :: S.IntArrayE
  , successes :: S.IntArrayE
  , predictors :: S.MatrixE
  , mDensity :: Maybe S.VectorE
  }

data BasicParameters = BasicParameters { mAlpha0 :: Maybe S.RealE
                                       , alpha :: S.VectorE
                                       , beta :: S.VectorE
                                       , logitMu :: S.VectorE
--                                       , mBetaDensity :: Maybe S.RealE
                                       }

modelData :: forall rs . (Typeable rs, F.ElemOf rs DT.PWPopPerSqMile)
          => S.DesignMatrixRow (F.Record rs)
          -> ModelConfig ()
          -> ACSRowTag rs (VU.Vector Int)
          -> S.StanModelBuilderEff [(F.Record rs, VU.Vector Int)] () (ModelData rs)
modelData dmr mc acsData = do
--  acsData <- S.dataSetTag @(F.Record rs, VU.Vector Int) S.ModelData "ACS"
  let nData' = S.dataSetSizeE acsData
      nStates' = S.groupSizeE stateGroup
  let trialsF v = v VU.! 0 + v VU.! 1
      successesF v = v VU.! 1
  trials' <- S.addCountData (modelIDT @rs @(VU.Vector Int)) acsData "trials" (trialsF . snd)
  successes' <- S.addCountData (modelIDT @rs @(VU.Vector Int)) acsData "successes" (successesF . snd)

  acsMat' <- S.addDesignMatrix (modelIDT @rs @(VU.Vector Int)) acsData (contramap fst dmr) Nothing
  let (_, nPredictors') = S.designMatrixColDimBinding dmr Nothing
  mDensity' <- case includeDensity mc of
    False -> pure Nothing
    True -> do
      rawDensity <- S.addRealData (modelIDT @rs @(VU.Vector Int)) acsData "rawLogDensity" Nothing Nothing (DDP.safeLog . F.rgetField @DT.PWPopPerSqMile . fst)
      stdDensity <- S.inBlock S.SBTransformedData $ S.addFromCodeWriter $ do
        let m = S.mean rawDensity
            sd = S.sqrt (S.variance rawDensity)
        S.declareRHSNW (S.NamedDeclSpec "stdLogDensity" $ S.vectorSpec nData')
          $ (rawDensity `S.minusE` m) `S.divideE` sd
      pure $ Just stdDensity

  pure $ ModelData acsData nData' nStates' nPredictors' trials' successes' acsMat' mDensity'

basicParameters :: ModelConfig () -> ModelData rs -> S.StanModelBuilderEff [Row rs a] () BasicParameters
basicParameters mc md = do
  mAlpha0P <- case mc.includeAlpha0  of
    True -> Just
            <$> S.simpleParameterWA
            (S.NamedDeclSpec "alpha0" $ S.realSpec)
            (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))
    False -> pure Nothing

  sigmaAlphaP <- S.simpleParameterWA
             (S.NamedDeclSpec "sigmaAlpha" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)
             (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  alphaP <- case mc.alphaType of
    HCentered -> S.addCenteredHierarchical
                 (S.NamedDeclSpec "alpha" $ S.vectorSpec md.nStates)
                 (S.given (S.realE 0) :> sigmaAlphaP :> TNil)
                 S.normalS
    HNonCentered -> S.simpleNonCentered
                    (S.NamedDeclSpec "alpha" $ S.vectorSpec md.nStates)
                    S.TransformedParametersBlock
                    (S.vectorSpec md.nStates)
                    (S.DensityWithArgs S.normalS $ S.realE 0 :> S.realE 1 :> TNil)
                    (S.given (S.realE 0) :> sigmaAlphaP :> TNil)
                    (\(ma :> sa :> TNil) r -> ma `S.plusE` (sa `S.timesE` r))

  betaP <- S.simpleParameterWA
         (S.NamedDeclSpec "beta" $ S.vectorSpec md.nPredictors)
         (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 2 :> TNil))

  mBetaDensityP <- case includeDensity mc of
    True -> Just <$> S.simpleParameterWA
                   (S.NamedDeclSpec "beta_Density" $ S.realSpec)
                   (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))
    False -> pure Nothing

  let f = S.parameterExpr
--      by v i = S.indexE S.s0 i v
      mBetaDensityTerm = S.timesE <$> (f <$> mBetaDensityP) <*> md.mDensity
      logitMu' = let x = (f alphaP `S.by` (S.dataByGroupIndexE md.acsDataTag stateGroup))
                       |+| (md.predictors |*| f betaP)
                 in addTermMaybe (f <$> mAlpha0P) (\a e -> a |+| e)
                    $ addTermMaybe mBetaDensityTerm (\bd e -> bd |+| e) x
  pure $ BasicParameters (f <$> mAlpha0P) (f alphaP) (f betaP) logitMu'


normalModel :: forall rs . (Typeable rs, F.ElemOf rs DT.PWPopPerSqMile)
            => S.DesignMatrixRow (F.Record rs)
            -> ModelConfig ()
            -> ACSRowTag rs (VU.Vector Int)
            -> S.StanModelBuilderEff [Row rs (VU.Vector Int)] () ()
normalModel dmr mc acsTag = do
  -- data
  md <- modelData dmr mc acsTag

  -- transformed data
  (obsP, binomialSigma2) <- S.inBlock S.SBTransformedData $ S.addFromCodeWriter $ do
    oP <- S.declareRHSNW (S.NamedDeclSpec "obsP" $ S.vectorSpec md.nData)
          $ S.to_vector md.successes |./| S.to_vector md.trials
    bS <- S.declareRHSNW (S.NamedDeclSpec "binomialSigma" $ S.vectorSpec md.nData)
          $ oP |.*| (S.realE 1 |-| oP) |./| S.to_vector md.trials
    pure (oP, bS)

  -- parameters & priors
  bParams <- basicParameters mc md

  sigmaP <- S.simpleParameterWA
         (S.NamedDeclSpec "sigma" $ S.addVMs (S.Modifiers [S.lowerM $ S.realE 0]) S.realSpec)
         (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  let sigma0 = S.parameterExpr sigmaP
      mu = S.inv_logit bParams.logitMu
      sigma = S.sqrt $ binomialSigma2 |+| (sigma0 |*| sigma0)
      ps = mu :> sigma :> TNil

  -- model
  S.inBlock S.SBModel $ S.addFromCodeWriter $ S.addStmt $ S.sample obsP S.normal ps

  -- generated quantities
  let vSpec = S.vectorSpec md.nData
      tempPs = do
        mu' <- S.declareRHSNW (S.NamedDeclSpec "muV" vSpec) mu
        s <- S.declareRHSNW (S.NamedDeclSpec "sigmaV" vSpec) sigma
        return (mu', s)
      tempP = S.declareRHSNW (S.NamedDeclSpec "pV" vSpec) obsP

  let at x n = S.sliceE S.s0 n x
  S.generateLogLikelihood
    md.acsDataTag
    S.normalDist
    ((\(e, sig) n -> e `at` n :> sig `at` n :> TNil) <$> tempPs)
    ((\o n ->  o `at` n) <$> tempP)

  _ <- S.inBlock S.SBGeneratedQuantities $ S.splitToGroupVars dmr bParams.beta (Just "beta")
  _ <- S.generatePosteriorPrediction'
    md.acsDataTag
    (S.NamedDeclSpec "pObserved" $ S.array1Spec md.nData S.realSpec)
    (\f n -> S.familyRNG S.normalDist (f n))
    ((\(e, sig) n -> e `at` n :> sig `at` n :> TNil) <$> tempPs)
    (\n p -> md.trials `at` n `S.timesE` p)
  pure ()


betaBinomialModel :: forall rs. (Typeable rs, F.ElemOf rs DT.PWPopPerSqMile)
                  => S.DesignMatrixRow (F.Record rs)
                  -> ModelConfig ()
                  -> ACSRowTag rs (VU.Vector Int)
                  -> S.StanModelBuilderEff [Row rs (VU.Vector Int)] () ()
betaBinomialModel dmr mc acsTag = do
  md <- modelData dmr mc acsTag
  absPredictors <- S.inBlock S.SBTransformedData $ S.addFromCodeWriter
                   $ S.declareRHSNW (S.NamedDeclSpec "absACSMat" $ S.matrixSpec md.nData md.nPredictors)
                   $ S.abs md.predictors
  -- parameters
  bParams <- basicParameters mc md
  let at x n = S.sliceE S.s0 n x
--      by v i = S.indexE S.s0 i v
  phiP <- S.addTransformedHP
          (S.NamedDeclSpec "phi" $ S.vectorSpec md.nPredictors)
          S.TransformedParametersBlock
          (Just $ S.Modifiers [S.lowerM $ S.realE 0, S.upperM $ S.realE 1]) -- constraints on phi_raw
          (S.DensityWithArgs S.betaS (S.realE 99 :> S.realE 1 :> TNil)) -- phi_raw is beta distributed
          (\t -> t |./| (S.realE 1 `S.minusE` t)) -- phi = phi_raw / (1 - phi_raw), component-wise

  let phi = S.parameterExpr phiP
      vSpec = S.vectorSpec md.nData
      tempPs = do
        mu <- S.declareRHSNW (S.NamedDeclSpec "muV" vSpec) $ S.inv_logit bParams.logitMu
        phiV <- S.declareRHSNW (S.NamedDeclSpec "mV" vSpec) $  absPredictors |*| phi
        betaA <- S.declareRHSNW (S.NamedDeclSpec "aV" vSpec) $ phiV |.*| mu
        betaB <-S.declareRHSNW (S.NamedDeclSpec "bV" vSpec) $ phiV |.*| (S.realE 1 |-| mu)
        pure (betaA, betaB)

  S.inBlock S.SBModel $ S.addFromCodeWriter $ do
    (betaA, betaB) <- tempPs
    let ps = md.trials :> betaA :> betaB :> TNil
    S.addStmt $ S.target $ S.densityE S.beta_binomial_lpmf md.successes ps

  S.generateLogLikelihood
    md.acsDataTag
    (S.betaBinomialDist' True)
    ((\(a, b) n -> md.trials `at` n :> a `at` n :> b `at` n :> TNil) <$> tempPs)
    (pure $ (md.successes `at`))

  _ <- S.inBlock S.SBGeneratedQuantities $ S.splitToGroupVars dmr bParams.beta (Just "beta")
  _ <- S.inBlock S.SBGeneratedQuantities $ S.splitToGroupVars dmr phi (Just "phi")
  _ <- S.generatePosteriorPrediction
    md.acsDataTag
    (S.NamedDeclSpec "pObserved" $ S.array1Spec md.nData S.intSpec)
    (S.betaBinomialDist' True)
    ((\(a, b) n -> md.trials `at` n :> a `at` n :> b `at` n :> TNil) <$> tempPs)
  pure ()

groupBuilderState :: forall rs a . (F.ElemOf rs GT.StateAbbreviation, Typeable rs, Typeable a)
                  => [Text]
                  -> S.StanDataBuilderEff S.ModelDataT [(F.Record rs, a)] (ACSRowTag rs a)
groupBuilderState states = do
  acsData <- S.addData "ACS" (modelIDT @rs @a) (S.ToFoldable id)
  S.addGroupIndexForData (modelIDT @rs @a) stateGroup acsData $ S.makeIndexFromFoldable show (F.rgetField @GT.StateAbbreviation . fst) states
  S.addGroupIntMapForData (modelIDT @rs @a) stateGroup acsData $ S.dataToIntMapFromFoldable (F.rgetField @GT.StateAbbreviation . fst) states
  pure acsData

{-
groupBuilderCD :: [Text] -> [Text] -> S.StanGroupBuilderM (F.FrameRec DDP.ACSByCD) () ()
groupBuilderCD states cds = do
  acsData <- S.addModelDataToGroupBuilder "ACS" (S.ToFoldable id)
  S.addGroupIndexForData stateGroup acsData $ S.makeIndexFromFoldable show (F.rgetField @GT.StateAbbreviation) states
  S.addGroupIndexForData cdGroup acsData $ S.makeIndexFromFoldable show DDP.districtKey cds
-}
cdGroup :: S.GroupTypeTag Text
cdGroup = S.GroupTypeTag "CD"

stateGroup :: S.GroupTypeTag Text
stateGroup = S.GroupTypeTag "State"

dmrS_ER :: forall rs . (F.ElemOf rs DT.Education4C
                       , F.ElemOf rs DT.SexC
                       , F.ElemOf rs DT.Race5C
                       )
                        => S.DesignMatrixRow (F.Record rs)
dmrS_ER = S.DesignMatrixRow "S_ER" [sexRP, raceEduRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    raceEduRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.R5_WhiteNonHispanic, DT.E4_HSGrad)) "RaceEdu"
                $ \r -> S.BEProduct2 (F.rgetField @DT.Race5C  r, F.rgetField @DT.Education4C r)

dmrC_S_ER :: forall rs . (F.ElemOf rs DT.CitizenC
                                   , F.ElemOf rs DT.Education4C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.Race5C
                                   )
                   => S.DesignMatrixRow (F.Record rs)
dmrC_S_ER = S.DesignMatrixRow "C_S_ER" [citRP, sexRP, raceEduRP]
  where
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    raceEduRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.R5_WhiteNonHispanic, DT.E4_HSGrad)) "RaceEdu"
                $ \r -> S.BEProduct2 (F.rgetField @DT.Race5C  r, F.rgetField @DT.Education4C r)

dmrS_CR :: forall rs . (F.ElemOf rs DT.CitizenC
                           , F.ElemOf rs DT.SexC
                           , F.ElemOf rs DT.Race5C
                           )
                   => S.DesignMatrixRow (F.Record rs)
dmrS_CR = S.DesignMatrixRow "S_CR" [sexRP, citRaceRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    citRaceRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.Citizen,  DT.R5_WhiteNonHispanic)) "CttRace"
                $ \r -> S.BEProduct2 (r ^. DT.citizenC, r ^. DT.race5C)

dmrC_S_A2R :: forall rs . (F.ElemOf rs DT.CitizenC
                         , F.ElemOf rs DT.SimpleAgeC
                         , F.ElemOf rs DT.SexC
                         , F.ElemOf rs DT.Race5C
                         )
                   => S.DesignMatrixRow (F.Record rs)
dmrC_S_A2R = S.DesignMatrixRow "C_S_A2R" [citRP, sexRP, ageRaceRP]
  where
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    ageRaceRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.Under, DT.R5_WhiteNonHispanic)) "Age2Race"
                $ \r -> S.BEProduct2 (r ^. DT.simpleAgeC, r ^. DT.race5C)


dmrS_A2ER :: forall rs . (F.ElemOf rs DT.Education4C
                            , F.ElemOf rs DT.SexC
                            , F.ElemOf rs DT.Race5C
                            , F.ElemOf rs DT.SimpleAgeC
                            )
                        => S.DesignMatrixRow (F.Record rs)
dmrS_A2ER = S.DesignMatrixRow "S_A2ER" [sexRP, ageRaceEduRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    ageRaceEduRP = S.boundedEnumRowPart (Just $ S.BEProduct3 (DT.Under, DT.R5_WhiteNonHispanic, DT.E4_HSGrad)) "Age2RaceEdu"
                   $ \r -> S.BEProduct3 (r ^. DT.simpleAgeC, r ^. DT.race5C, r ^. DT.education4C)


dmrS_A2CR :: forall rs . (F.ElemOf rs DT.CitizenC
                        , F.ElemOf rs DT.SexC
                        , F.ElemOf rs DT.Race5C
                        , F.ElemOf rs DT.Age5C
                        )
              => S.DesignMatrixRow (F.Record rs)
dmrS_A2CR = S.DesignMatrixRow "S_A2CR" [sexRP, citAgeRaceRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    citAgeRaceRP = S.boundedEnumRowPart (Just $ S.BEProduct3 (DT.Citizen, DT.Under, DT.R5_WhiteNonHispanic)) "CitAgeRace"
                   $ \r -> S.BEProduct3 (r ^. DT.citizenC, DT.age5ToSimple $ r ^. DT.age5C, r ^. DT.race5C)


dmrS_AR :: forall rs . (F.ElemOf rs DT.SexC
                       , F.ElemOf rs DT.Race5C
                       , F.ElemOf rs DT.Age5C
                       )
              => S.DesignMatrixRow (F.Record rs)
dmrS_AR = S.DesignMatrixRow "S_AR" [sexRP, ageRaceRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    ageRaceRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.A5_35To44, DT.R5_WhiteNonHispanic)) "AgeRace"
                $ \r -> S.BEProduct2 (r ^. DT.age5C, r ^. DT.race5C)


dmrS_CAR :: forall rs . (F.ElemOf rs DT.CitizenC
                        , F.ElemOf rs DT.SexC
                        , F.ElemOf rs DT.Race5C
                        , F.ElemOf rs DT.Age5C
                        )
              => S.DesignMatrixRow (F.Record rs)
dmrS_CAR = S.DesignMatrixRow "S_CAR" [sexRP, citAgeRaceRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    citAgeRaceRP = S.boundedEnumRowPart (Just $ S.BEProduct3 (DT.Citizen, DT.A5_35To44, DT.R5_WhiteNonHispanic)) "CitAgeRace"
                   $ \r -> S.BEProduct3 (r ^. DT.citizenC, r ^. DT.age5C, r ^. DT.race5C)


dmrS_C_AR :: forall rs . (F.ElemOf rs DT.CitizenC
                        , F.ElemOf rs DT.SexC
                        , F.ElemOf rs DT.Race5C
                        , F.ElemOf rs DT.Age5C
                        )
              => S.DesignMatrixRow (F.Record rs)
dmrS_C_AR = S.DesignMatrixRow "S_C_AR" [sexRP, citRP, ageRaceRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    ageRaceRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.A5_35To44, DT.R5_WhiteNonHispanic)) "AgeRace"
                $ \r -> S.BEProduct2 (r ^. DT.age5C, r ^. DT.race5C)


designMatrixRowAge :: forall rs . (F.ElemOf rs DT.CitizenC
                                  , F.ElemOf rs DT.Education4C
                                  , F.ElemOf rs DT.SexC
                                  , F.ElemOf rs DT.Race5C
                                  )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowAge = S.DesignMatrixRow "DMAge" [citRP, sexRP, eduRP, raceRP]
  where
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    eduRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Education" (F.rgetField @DT.Education4C)
    raceRP = S.boundedEnumRowPart (Just DT.R5_WhiteNonHispanic) "Race" (F.rgetField @DT.Race5C)



designMatrixRowCitizen :: forall rs . (F.ElemOf rs DT.Education4C
                                      , F.ElemOf rs DT.SexC
                                      , F.ElemOf rs DT.Race5C
                                      )
                       => S.DesignMatrixRow (F.Record rs)
designMatrixRowCitizen = S.DesignMatrixRow "DMCitizen" [sexRP, eduRP, raceRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    eduRP = S.boundedEnumRowPart (Just DT.E4_HSGrad) "Education" (F.rgetField @DT.Education4C)
    raceRP = S.boundedEnumRowPart (Just DT.R5_WhiteNonHispanic) "Race" (F.rgetField @DT.Race5C)



--
newtype ModelResult2 ks = ModelResult2 { unModelResult :: Map (F.Record ks) Double }

deriving stock instance (Show (F.Record ks)) => Show (ModelResult2 ks)

instance (V.RMap ks, FS.RecFlat ks, Ord (F.Rec FS.SElField ks), Ord (F.Record ks)) => Flat.Flat (ModelResult2 ks) where
  size = Flat.size . M.mapKeys FS.toS . unModelResult
  encode = Flat.encode . M.mapKeys FS.toS . unModelResult
  decode = fmap (ModelResult2 . M.mapKeys FS.fromS) Flat.decode


applyModelResult2 :: (ks F.⊆ rs, Ord (F.Record ks), Show (F.Record rs))
                 => ModelResult2 ks -> F.Record rs -> Either Text Double
applyModelResult2 (ModelResult2 m) r = case M.lookup (F.rcast r) m of
                                         Nothing -> Left $ "applyModelResult2: key=" <> show r <> " not found in model result map."
                                         Just p -> Right p



--
data ModelResult g ks = ModelResult { alpha0 :: Double, geoAlpha :: Map g Double, ldSI :: (Double, Double), catAlpha :: Map (F.Record ks) Double }

deriving stock instance (Show (F.Record ks), Show g) => Show (ModelResult g ks)

modelResultToFTuple :: (Ord (V.Rec FS.SElField ks), V.RMap ks) => ModelResult g ks -> (Double, Map g Double, (Double, Double), Map (V.Rec FS.SElField ks) Double)
modelResultToFTuple (ModelResult a b c d) = (a, b, c, M.mapKeys FS.toS d)

modelResultFromFTuple :: (Ord (F.Record ks), V.RMap ks) => (Double, Map g Double, (Double, Double), Map (V.Rec FS.SElField ks) Double) -> ModelResult g ks
modelResultFromFTuple (a, b, c, d) = ModelResult a b c (M.mapKeys FS.fromS d)

instance (V.RMap ks, FS.RecFlat ks, Flat.Flat g, Ord g, Ord (F.Rec FS.SElField ks), Ord (F.Record ks)) => Flat.Flat (ModelResult g ks) where
  size = Flat.size . modelResultToFTuple
  encode = Flat.encode . modelResultToFTuple
  decode = fmap modelResultFromFTuple Flat.decode


applyModelResult :: (F.ElemOf rs DT.PWPopPerSqMile, ks F.⊆ rs, Ord g, Show g, Ord (F.Record ks), Show (F.Record rs))
                 => ModelResult g ks -> g -> F.Record rs -> Either Text Double
applyModelResult (ModelResult a ga (ldS, ldI) ca) g r = invLogit <$> xE where
  invLogit y = 1 / (1 + Numeric.exp (negate y))
  geoXE = maybe (Left $ "applyModelResult: " <> show g <> " missing from geography alpha map") Right $ M.lookup g ga
  densX = ldI + ldS * (DDP.safeLog $ F.rgetField @DT.PWPopPerSqMile r)
  catXE = maybe (Left $ "applyModelResult: " <> show r <> " missing from category alpha map") Right $ M.lookup (F.rcast r) ca
  xE = (\a' d g' c -> a' + d + g' + c) <$> pure a <*> pure densX <*> geoXE <*> catXE

stateModelResultAction :: forall rs ks a r gq.
                          (K.KnitEffects r
                          , Typeable rs
                          , Typeable a
                          , F.ElemOf rs DT.PWPopPerSqMile
--                          , ks F.⊆ rs
                          , Ord (F.Record ks)
                          , BRK.FiniteSet (F.Record ks)
                          )
                       => ModelConfig Text
                       -> S.DesignMatrixRow (F.Record ks)
                       -> S.ResultAction [(F.Record rs, a)] gq S.DataSetGroupIntMaps S.DataSetGroupIntMaps r () (ModelResult Text ks)
stateModelResultAction mc dmr = S.UseSummary f where
  f summary _ modelDataAndIndexes_C _ = do
--    let resultCacheKey = modelID mc <> "_" <> S.dmName dmr <> modelConfigSuffix mc
    (modelData', resultIndexesE) <- K.ignoreCacheTime modelDataAndIndexes_C
    -- we need to rescale the density component to work
    let premap = DDP.safeLog . F.rgetField @DT.PWPopPerSqMile . fst
        msFld = (,) <$> FL.mean <*> FL.std
        (ldMean, ldSigma) = FL.fold (FL.premap premap msFld) modelData'
    stateIM <- K.knitEither
      $ resultIndexesE >>= S.getGroupIndex (S.RowTypeTag @(Row rs a) "ACS") stateGroup
    let getScalar n = K.knitEither $ S.getScalar . fmap CS.mean <$> S.parseScalar n (CS.paramStats summary)
        getVector n = K.knitEither $ S.getVector . fmap CS.mean <$> S.parse1D n (CS.paramStats summary)
    alpha' <- case mc.includeAlpha0  of
      False -> pure 0
      True -> getScalar "alpha0"
    geoMap <- (\stIM alphaV -> M.fromList $ zip (IM.elems stIM) (Vec.toList alphaV)) <$> pure stateIM <*> getVector "alpha"
    (ldSlope, ldIntercept) <- case includeDensity mc of
      False -> pure (0, 0)
      True -> (\x -> (x / ldSigma, negate $ x * ldMean / ldSigma)) <$> getScalar "beta_Density"
    catBeta <- VU.convert <$> getVector "beta"
    let (S.MatrixRowFromData _ _ _ rowVecF) = S.matrixFromRowData dmr Nothing
        allCatRows = Set.toList $ BRK.elements @(F.Record ks)
        g v1 v2 = VU.foldl' (\a (b, c) -> a + b * c) 0 $ VU.zip v1 v2
        catMap = M.fromList $ zip allCatRows (g catBeta . rowVecF <$> allCatRows)
    pure $ ModelResult alpha' geoMap (ldSlope, ldIntercept) catMap

--    modelResult <- ModelResult <$> getVector "alpha"
type CitizenStateModelResult = ModelResult Text [DT.SexC, DT.Education4C, DT.Race5C]
type AgeStateModelResult = ModelResult Text [DT.CitizenC, DT.SexC, DT.Education4C, DT.Race5C]
type EduStateModelResult = ModelResult Text [DT.Age5C, DT.SexC, DT.Race5C]


logDensityDMRP :: F.ElemOf rs DT.PWPopPerSqMile => S.DesignMatrixRowPart (F.Record rs)
logDensityDMRP = S.DesignMatrixRowPart "Density" 1 DDP.logDensityPredictor

----

categoricalModel :: forall rs . Typeable rs
                 => Int
                 -> S.DesignMatrixRow (F.Record rs)
                 -> ACSRowTag rs (VU.Vector Int)
                 -> S.StanModelBuilderEff [Row rs (VU.Vector Int)] () ()
categoricalModel numInCat dmr acsData = do
--  acsData <- S.dataSetTag @(F.Record rs, VU.Vector Int) S.ModelData "ACS"
  let nDataE = S.dataSetSizeE acsData
  nInCatE <- S.addFixedIntModel @[Row rs (VU.Vector Int)] "K" numInCat
  countsE <- S.addIntArrayData (modelIDT @rs @(VU.Vector Int)) acsData "counts" nInCatE (Just 0) Nothing snd
  acsMatE <- S.addDesignMatrix (modelIDT @rs @(VU.Vector Int)) acsData (contramap fst dmr) Nothing
  let (_, nPredictorsE) = S.designMatrixColDimBinding dmr Nothing
  -- parameters
  -- zero vector for identifiability trick
  zvP <- S.addBuildParameter
         $ S.TransformedDataP
         $ S.TData
         (S.NamedDeclSpec "zeroes" $ S.vectorSpec nPredictorsE)
         []
         TNil
         (const $ S.DeclRHS $ S.rep_vector (S.realE 0) nPredictorsE)

  betaRawP <- S.addBuildParameter
              $ S.UntransformedP
              (S.NamedDeclSpec "beta_raw" $ S.matrixSpec nPredictorsE (nInCatE |-| S.intE 1))
              []
              TNil
              (\_ _ -> pure ())

  betaP <- S.addBuildParameter
           $ S.TransformedP
           (S.NamedDeclSpec "beta" $ S.matrixSpec nPredictorsE nInCatE)
           []
           (betaRawP :> zvP :> TNil)
           S.TransformedParametersBlock
           (\(betaRawE :> zvE :> TNil) -> S.DeclRHS $ S.append_col betaRawE zvE)
           (S.given (S.realE 0) :> S.given (S.realE 2) :> TNil)
           (\normalPS x -> S.addStmt $ S.sample (S.to_vector x) S.normalS normalPS)

  let betaE = S.parameterExpr betaP
      betaXD = S.declareRHSNW
               (S.NamedDeclSpec "beta_x" $ S.matrixSpec nDataE nInCatE)
               (acsMatE `S.timesE` betaE)
      at x n = S.sliceE S.s0 n x

  S.inBlock S.SBModel $ S.addStmtToCode $ S.cwStmt_ $ do
--    let sizeE e = S.functionE S.size (e :> TNil)
    betaX <- betaXD
    S.addStmt $ S.for "n" (S.SpecificNumbered (S.intE 1) nDataE) $ \n ->
      S.target $ S.densityE S.multinomial_logit_lupmf (countsE `at` n) (S.transposeE (betaX `at` n) :> TNil)

  gqBetaX <- S.inBlock S.SBLogLikelihood $ S.addFromCodeWriter betaXD
  S.generateLogLikelihood
    acsData
    S.multinomialLogitDist
    (pure $ \nE -> S.transposeE (gqBetaX `at` nE) :> TNil)
    (pure $ \nE -> countsE `at` nE)


designMatrixRowEdu3 :: forall rs . (F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.RaceAlone4C
                                   , F.ElemOf rs DT.HispC
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu3 = S.DesignMatrixRow "DMEdu3" [sexRaceAgeRP]
  where
--    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC . fst)
    race5Census r = DT.race5FromRaceAlone4AndHisp True (F.rgetField @DT.RaceAlone4C r) (F.rgetField @DT.HispC r)
    sexRaceAgeRP = S.boundedEnumRowPart (Just $ S.BEProduct3 (DT.Female, DT.R5_WhiteNonHispanic, DT.A5_35To44)) "SexRaceAge"
                $ \r -> S.BEProduct3 (F.rgetField @DT.SexC r, race5Census r, F.rgetField @DT.Age5C r)

designMatrixRowEdu4 :: forall rs . (F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.RaceAlone4C
                                   , F.ElemOf rs DT.HispC
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu4 = S.DesignMatrixRow "DMEdu4" [sexRaceAgeRP]
  where
    race5Census r = DT.race5FromRaceAlone4AndHisp True (F.rgetField @DT.RaceAlone4C r) (F.rgetField @DT.HispC r)
    sexRaceAgeRP = S.boundedEnumRowPart Nothing "SexRaceAge"
                $ \r -> S.BEProduct3 (F.rgetField @DT.SexC r, race5Census r, F.rgetField @DT.Age5C r)

designMatrixRowEdu8 :: forall rs . (F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.RaceAlone4C
                                   , F.ElemOf rs DT.HispC
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu8 = S.DesignMatrixRow "DMEdu8" [sexRP, raceAgeRP]
  where
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    race5Census r = DT.race5FromRaceAlone4AndHisp True (F.rgetField @DT.RaceAlone4C r) (F.rgetField @DT.HispC r)
    raceAgeRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.R5_WhiteNonHispanic, DT.A5_35To44)) "RaceAge"
                $ \r -> S.BEProduct2 (race5Census r, F.rgetField @DT.Age5C r)

designMatrixRowEdu5 :: forall rs . ( F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.RaceAlone4C
                                   , F.ElemOf rs DT.HispC
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu5 = S.DesignMatrixRow "DMEdu5" [sexRP, ageRP, raceRP, sexRaceAgeRP]
  where
    race5Census r = DT.race5FromRaceAlone4AndHisp True (F.rgetField @DT.RaceAlone4C r) (F.rgetField @DT.HispC r)
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    ageRP = S.boundedEnumRowPart (Just  DT.A5_35To44) "Age" (F.rgetField @DT.Age5C)
    raceRP = S.boundedEnumRowPart (Just DT.R5_WhiteNonHispanic) "Race" race5Census
    sexRaceAgeRP = S.boundedEnumRowPart (Just $ S.BEProduct3 (DT.Female, DT.R5_WhiteNonHispanic, DT.A5_35To44)) "SexRaceAge"
                $ \r -> S.BEProduct3 (F.rgetField @DT.SexC r, race5Census r, F.rgetField @DT.Age5C r)

designMatrixRowEdu6 :: forall rs . (F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.RaceAlone4C
                                   , F.ElemOf rs DT.HispC
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu6 = S.DesignMatrixRow "DMEdu6" [sexRP, ageRP, raceRP, sexRaceAgeRP]
  where
    race5Census r = DT.race5FromRaceAlone4AndHisp True (F.rgetField @DT.RaceAlone4C r) (F.rgetField @DT.HispC r)
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    ageRP = S.boundedEnumRowPart Nothing "Age" (F.rgetField @DT.Age5C)
    raceRP = S.boundedEnumRowPart Nothing "Race" race5Census
    sexRaceAgeRP = S.boundedEnumRowPart Nothing "SexRaceAge"
                $ \r -> S.BEProduct3 (F.rgetField @DT.SexC r, race5Census r, F.rgetField @DT.Age5C r)


designMatrixRowEdu :: forall rs . (F.ElemOf rs DT.CitizenC
                                  ,  F.ElemOf rs DT.Age5C
                                  ,  F.ElemOf rs DT.SexC
                                  ,  F.ElemOf rs DT.Race5C
                                  )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu = S.DesignMatrixRow "DMEdu" [citRP, sexRP, ageRP, raceRP]
  where
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC )
    ageRP = S.boundedEnumRowPart (Just DT.A5_35To44) "Age" (F.rgetField @DT.Age5C)
    raceRP = S.boundedEnumRowPart (Just DT.R5_WhiteNonHispanic) "Race" (F.rgetField @DT.Race5C)

designMatrixRowEdu7 :: forall rs . (F.ElemOf rs DT.CitizenC
                                   , F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.Race5C
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu7 = S.DesignMatrixRow "DMEdu7" [citRP, sexRP, raceAgeRP]
  where
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    raceAgeRP = S.boundedEnumRowPart Nothing "RaceAge"
                $ \r -> S.BEProduct2 (F.rgetField @DT.Race5C r, F.rgetField @DT.Age5C r)

designMatrixRowEdu2 :: forall rs . (F.ElemOf rs DT.CitizenC
                                   , F.ElemOf rs DT.Age5C
                                   , F.ElemOf rs DT.SexC
                                   , F.ElemOf rs DT.Race5C
                                   )
                   => S.DesignMatrixRow (F.Record rs)
designMatrixRowEdu2 = S.DesignMatrixRow "DMEdu2" [citRP, sexRP, raceAgeRP]
  where
    citRP = S.boundedEnumRowPart Nothing "Citizen" (F.rgetField @DT.CitizenC )
    sexRP = S.boundedEnumRowPart Nothing "Sex" (F.rgetField @DT.SexC)
    raceAgeRP = S.boundedEnumRowPart (Just $ S.BEProduct2 (DT.R5_WhiteNonHispanic, DT.A5_35To44)) "RaceAge"
                $ \r -> S.BEProduct2 (F.rgetField @DT.Race5C r, F.rgetField @DT.Age5C r)




{-
binomialNormalModel :: forall rs . Typeable rs
                 => S.DesignMatrixRow (F.Record rs)
                 -> S.StanBuilderM [(F.Record rs, VU.Vector Int)] () ()
binomialNormalModel dmr = do
  acsData <- S.dataSetTag @(F.Record rs, VU.Vector Int) S.ModelData "ACS"
  let nData = S.dataSetSizeE acsData
      nStates = S.groupSizeE stateGroup
--  countsE <- S.addIntArrayData acsData "counts" (S.intE 2) (Just 0) Nothing snd
  let trials v = v VU.! 0 + v VU.! 1
      successes v = v VU.! 1
  trials <- S.addCountData acsData "trials" (trials . snd)
  successes <- S.addCountData acsData "successes" (successes . snd)
  acsMat <- S.addDesignMatrix acsData (contramap fst dmr) Nothing
  let (_, nPredictors) = S.designMatrixColDimBinding dmr Nothing
      at x n = S.sliceE S.s0 n x
  -- parameters
  sigmaAlphaP <- S.simpleParameterWA
             (S.NamedDeclSpec "sigmaAlpha" $ S.realSpec [S.lowerM $ S.realE 0])
             (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  alphaP <- S.addCenteredHierarchical
            (S.NamedDeclSpec "alpha" $ S.vectorSpec nStates [])
            (S.given (S.realE 0) :> S.build sigmaAlphaP :> TNil)
            S.normalS

  betaP <- S.simpleParameterWA
           (S.NamedDeclSpec "beta" $ S.vectorSpec nPredictors [])
           (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 2 :> TNil))

  muErrP <- S.simpleParameterWA
            (S.NamedDeclSpec "muErr" $ S.realSpec [])
            (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  sigmaErrP <- S.simpleParameterWA
               (S.NamedDeclSpec "sigmaErr" $ S.realSpec [S.lowerM $ S.realE 0])
               (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  errP <- S.addCenteredHierarchical
          (S.NamedDeclSpec "err" $ S.vectorSpec nData [])
          (S.build muErrP :> S.build sigmaErrP :> TNil)
          S.normalS

  let alpha = S.parameterTagExpr alphaP
      beta  = S.parameterTagExpr betaP
      err = S.parameterTagExpr errP
      p =  S.indexE S.s0 (S.byGroupIndexE acsData stateGroup) alpha `S.plusE` (acsMat `S.timesE` beta) `S.plusE` err
      vSpec = S.vectorSpec nData []
      tmpP = S.declareRHSNW (S.NamedDeclSpec "pV" vSpec) p

  S.inBlock S.SBModel $ S.addFromCodeWriter $ do
    p <- tmpP
    S.addStmt $ S.for "n" (S.SpecificNumbered (S.intE 1) nData) $ \n ->
      let lhs = successes `at` n
          ps = trials `at` n :> p `at` n :> TNil
      in [S.target $ S.densityE S.binomial_logit_lpmf lhs ps]

  S.generateLogLikelihood
    acsData
    S.binomialLogitDist
    ((\p n -> trials `at` n :> p `at` n :> TNil) <$> tmpP)
    (pure (successes `at`))

  _ <- S.generatePosteriorPrediction
    acsData
    (S.NamedDeclSpec "pObserved" $ S.array1Spec nData $ S.intSpec [])
    S.binomialLogitDist
    ((\p n -> trials `at` n :> p `at` n :> TNil) <$> tmpP)
  pure ()

binomialModel :: forall rs . Typeable rs
                 => S.DesignMatrixRow (F.Record rs)
                 -> S.StanBuilderM [(F.Record rs, VU.Vector Int)] () ()
binomialModel dmr = do
  acsData <- S.dataSetTag @(F.Record rs, VU.Vector Int) S.ModelData "ACS"
  let nData = S.dataSetSizeE acsData
      nStates = S.groupSizeE stateGroup
      trialsF v = v VU.! 0 + v VU.! 1
      successesF v = v VU.! 1
  trials <- S.addCountData acsData "trials" (trialsF . snd)
  successes <- S.addCountData acsData "successes" (successesF . snd)
  acsMat <- S.addDesignMatrix acsData (contramap fst dmr) Nothing
  let (_, nPredictors) = S.designMatrixColDimBinding dmr Nothing
      at x n = S.sliceE S.s0 n x
      by v i = S.indexE S.s0 i v

  -- parameters

  betaP <- S.simpleParameterWA
           (S.NamedDeclSpec "beta" $ S.vectorSpec nPredictors [])
           (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 2 :> TNil))

  sigmaAlphaP <- S.simpleParameterWA
             (S.NamedDeclSpec "sigmaAlpha" $ S.realSpec [S.lowerM $ S.realE 0])
             (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  alphaP <- S.addCenteredHierarchical
            (S.NamedDeclSpec "alpha" $ S.vectorSpec nStates [])
            (S.given (S.realE 0) :> S.build sigmaAlphaP :> TNil)
            S.normalS

  let beta = S.parameterTagExpr betaP
      alpha = S.parameterTagExpr alphaP
      logitMu = (alpha `by` (S.byGroupIndexE acsData stateGroup)) `S.plusE` (acsMat `S.timesE` beta)
      vSpec = S.vectorSpec nData []
      tempLM = S.declareRHSNW (S.NamedDeclSpec "lmV" vSpec) logitMu

  S.inBlock S.SBModel $ S.addFromCodeWriter $ do
    S.addStmt $ S.target $ S.densityE S.binomial_logit_lpmf successes (trials :> logitMu :> TNil)

  S.generateLogLikelihood
    acsData
    S.binomialLogitDist
    ((\lm n -> (trials `at` n :> lm `at` n :> TNil)) <$> tempLM)
    (pure $ \n -> successes `at` n)

  S.inBlock S.SBGeneratedQuantities $ S.splitToGroupVars dmr beta (Just "beta")
  _ <- S.generatePosteriorPrediction
    acsData
    (S.NamedDeclSpec "pObserved" $ S.array1Spec nData $ S.intSpec [])
    (S.binomialLogitDist' True)
    ((\lm n -> trials `at` n :> lm `at` n :> TNil) <$> tempLM)
  pure ()
-}


{-
negBinomialModel :: forall rs.Typeable rs
                 => S.DesignMatrixRow (F.Record rs)
                 -> S.StanBuilderM [(F.Record rs, VU.Vector Int)] () ()
negBinomialModel dmr = do
  acsData <- S.dataSetTag @(F.Record rs, VU.Vector Int) S.ModelData "ACS"
  let nData = S.dataSetSizeE acsData
      nStates = S.groupSizeE stateGroup
      trials v = v VU.! 0 + v VU.! 1
      successes v = v VU.! 1
  trials <- S.addCountData acsData "trials" (trials . snd)
  successes <- S.addCountData acsData "successes" (successes . snd)
  acsMat <- S.addDesignMatrix acsData (contramap fst dmr) Nothing
  let (_, nPredictors) = S.designMatrixColDimBinding dmr Nothing
      at x n = S.sliceE S.s0 n x
      by v d g = TE,indexE S.s0 (S.byGroupIndexE d g) v
      vSpec = S.vectorSpec nData []
  -- transformed data
  realSuccesses <- S.inBlock S.SBTransformedData $ S.addFromCodeWriter
                   $ S.declareRHSNW (S.NamedDeclSpec "rSuccesses" vSpec)
                   $ S.functionE S.to_vector (successes :> TNil)

  -- parameters
  sigmaAlphaP <- S.simpleParameterWA
             (S.NamedDeclSpec "sigmaAlpha" $ S.realSpec [S.lowerM $ S.realE 0])
             (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  alphaP <- S.addCenteredHierarchical
            (S.NamedDeclSpec "alpha" $ S.vectorSpec nStates [])
            (S.given (S.realE 0) :> S.build sigmaAlphaP :> TNil)
            S.normalS

  betaP <- S.simpleParameterWA
           (S.NamedDeclSpec "beta" $ S.vectorSpec nPredictors [])
           (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 2 :> TNil))

  phiP <- S.simpleParameterWA
           (S.NamedDeclSpec "phi" $ S.vectorSpec nPredictors [])
           (S.DensityWithArgs S.normalS (S.realE 0 :> S.realE 1 :> TNil))

  let alpha = S.parameterTagExpr alphaP
      beta  = S.parameterTagExpr betaP
      phi = S.parameterTagExpr phiP
      p =  (alpha `by` acsData stateGroup) `S.plusE` (acsMat `S.timesE` beta)
      tmpMu = S.declareRHSNW (S.NamedDeclSpec "muV" vSpec) $ p `eltTimes` realSuccesses

  S.inBlock S.SBModel $ S.addFromCodeWriter $ do
    mu <- tmpMu
    S.addStmt $ S.target $ S.densityE S.neg_binomial_2 trials (mu :> phi :> TNil)

  S.generateLogLikelihood
    acsData
    S.binomialLogitDist
    ((\p n -> trials `at` n :> p `at` n :> TNil) <$> tmpP)
    (pure (successes `at`))

  _ <- S.generatePosteriorPrediction
    acsData
    (S.NamedDeclSpec "pObserved" $ S.array1Spec nData $ S.intSpec [])
    S.binomialLogitDist
    ((\p n -> trials `at` n :> p `at` n :> TNil) <$> tmpP)
  pure ()
-}
