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

module BlueRipple.Model.Demographic.NullSpaceBasis
  (
    module BlueRipple.Model.Demographic.NullSpaceBasis
  )
where

import qualified BlueRipple.Model.Demographic.EnrichData as DED
import qualified Data.Map.Strict as M
import qualified Data.List as List

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Array as A
import qualified Numeric.LinearAlgebra.Array.Util as A
import qualified Data.Vector.Storable as VS

-- # of dimensions at each index
newtype Dimensions = Dimensions { dimensions :: [Int] }

-- # a set of indices into the array
newtype Indices = Indices { indices :: [Int] }

-- a subset specified by choosing some indices
newtype Subset = Subset { subset :: [Int] }

-- given a set of stencils, i.e., rows of constraints on distributions,
-- which are possibly linearly-dependent,
-- compute an orthonormal basis for them
-- as well as an orthonormal basis for their orthogonal complement, the null space
nullSpaceSVD :: Int -> [DED.Stencil Int] -> (LA.Matrix LA.R, LA.Matrix LA.R)
nullSpaceSVD n sts =
  let cM = DED.mMatrix n sts
      (u, s, v') = LA.svd cM
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
ds # cs = A.listArray ds cs :: A.Array Double

rArray :: Char -> LA.Vector LA.R -> Int -> A.Array LA.R
rArray i m p =
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
      rT k = rArray (indexLabels List.!! (k - 1)) (a List.!! (k - 1)) k
      rTs = fmap rT [2..nDims]
  in asStackedRowVec $ List.foldl' (A.|*|) (rT 1) rTs

interactionVC :: Dimensions -> Subset -> Indices -> LA.Vector LA.R
interactionVC dims sub subsetIs =
  let a = aVecs dims sub subsetIs
      nDims = length (dimensions dims)
      indexLabels = take nDims ['a'..]
      rT k = rArray (indexLabels List.!! (k - 1)) (a List.!! (k - 1)) k
      rTs = fmap rT [2..nDims]
  in asStackedColVec $ List.foldl' (A.|*|) (rT 1) rTs


subsetIndices :: Dimensions -> Subset -> [Indices]
subsetIndices dims subs = fmap Indices $ traverse (\si -> [1..((dimensions dims List.!! (si - 1)) - 1)]) $ subset subs

interactionBasis :: Dimensions -> Subset -> LA.Matrix LA.R
interactionBasis dims sub = LA.fromColumns $ fmap oneVec subIs
  where
    oneVec = interactionV dims sub
    subIs = subsetIndices dims sub
