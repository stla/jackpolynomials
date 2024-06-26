{-|
Module      : Math.Algebra.JackPol
Description : Symbolic Jack polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

Computation of symbolic Jack polynomials, zonal polynomials, Schur polynomials and skew Schur polynomials. 
See README for examples and references.
-}

{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.JackPol
  ( jackPol', skewJackPol', zonalPol', schurPol', skewSchurPol'
  , jackPol, skewJackPol, zonalPol, schurPol, skewSchurPol )
  where
import           Prelude 
  hiding ((*), (+), (-), (/), (^), (*>), product, sum, fromIntegral, fromInteger)
import           Algebra.Additive           ( (+), (-), sum )
import qualified Algebra.Field              as AlgField
import           Algebra.Ring               ( (*), product, one, fromInteger )
import qualified Algebra.Ring               as AlgRing
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import qualified Data.HashMap.Strict        as HM
import qualified Data.Map.Strict            as DM
import           Data.Maybe                 ( fromJust, isJust )
import           Math.Algebra.Jack.Internal ( _betaratio, jackCoeffC
                                            , _N, _isPartition, Partition
                                            , jackCoeffP, jackCoeffQ
                                            , skewSchurLRCoefficients
                                            , isSkewPartition, _fromInt
                                            , skewJackInMSPbasis )
import           Math.Algebra.Hspray        ( FunctionLike (..), (.^)
                                            , lone, Spray, QSpray
                                            , zeroSpray, unitSpray
                                            , fromList )
import           Math.Combinat.Permutations ( permuteMultiset )

-- | Symbolic Jack polynomial
jackPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Rational  -- ^ Jack parameter
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> QSpray
jackPol' = jackPol

-- | Symbolic Jack polynomial
jackPol :: forall a. (Eq a, AlgField.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> a         -- ^ Jack parameter
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> Spray a
jackPol n lambda alpha which 
  | n < 0 = error "jackPol: negative number of variables."
  | not (_isPartition lambda) = error "jackPol: invalid integer partition."
  | not (which `elem` ['J', 'C', 'P', 'Q']) = 
      error "jackPol: please use 'J', 'C', 'P' or 'Q' for last argument."
  | n == 0 = if null lambda
      then unitSpray
      else zeroSpray
  | otherwise =
    case which of 
      'J' -> resultJ
      'C' -> jackCoeffC lambda alpha *^ resultJ
      'P' -> jackCoeffP lambda alpha *^ resultJ
      _   -> jackCoeffQ lambda alpha *^ resultJ
      where
      resultJ = jac (length x) 0 lambda lambda arr0 one
      nll = _N lambda lambda
      x = map lone [1 .. n] :: [Spray a]
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then one
        else product [i .^ alpha + one | i <- [1 .. nu0-1]]
      jac :: Int -> Int -> Partition -> Partition 
              -> Array (Int,Int) (Maybe (Spray a)) -> a -> Spray a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = unitSpray
        | length nu > m && nu !! m > 0    = zeroSpray
        | m == 1                          = 
            theproduct (nu!!0) *^ (x!!0 ^**^ nu!!0) 
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (beta *^ (jac (m-1) 0 nu nu arr one ^*^ 
                      ((x!!(m-1)) ^**^ (sum mu - sum nu)))) (max 1 k)
            go :: Spray a -> Int -> Spray a
            go !ss ii
              | length nu < ii || nu!!(ii-1) == 0 = ss
              | otherwise =
                let u = nu!!(ii-1) in
                if length nu == ii && u > 0 || u > nu !! ii
                  then
                    let nu'   = (element (ii-1) .~ u-1) nu in
                    let gamma = beta * _betaratio mu nu ii alpha in
                    if u > 1
                      then
                        go (ss ^+^ jac m ii mu nu' arr gamma) (ii + 1)
                      else
                        if nu'!!0 == 0
                          then
                            go (ss ^+^ (gamma *^ (x!!(m-1) ^**^ sum mu))) 
                                (ii + 1)
                          else
                            let arr' = arr // [((_N lambda nu, m), Just ss)] in
                            let jck  = jac (m-1) 0 nu' nu' arr' one in
                            let jck' = gamma *^ (jck ^*^ 
                                        (x!!(m-1) ^**^ (sum mu - sum nu'))) in
                            go (ss ^+^ jck') (ii + 1)
                  else
                    go ss (ii + 1)

skewJackPol :: 
    (Eq a, AlgField.C a) 
  => Int 
  -> Partition 
  -> Partition 
  -> a
  -> Char 
  -> Spray a
skewJackPol n lambda mu alpha which = 
  HM.unions sprays
  where
    msCombo = 
      DM.filterWithKey 
        (\kappa _ -> length kappa <= n) 
          (skewJackInMSPbasis alpha which lambda mu)
    sprays = 
      map (
        \(kappa, coeff) -> 
          fromList
            (zip 
              (permuteMultiset (kappa ++ replicate (n - length kappa) 0)) 
              (repeat coeff))
        ) (DM.assocs msCombo)

skewJackPol' :: 
     Int 
  -> Partition 
  -> Partition 
  -> Rational
  -> Char 
  -> QSpray
skewJackPol' = skewJackPol

-- | Symbolic zonal polynomial
zonalPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> QSpray
zonalPol' = zonalPol

-- | Symbolic zonal polynomial
zonalPol :: (Eq a, AlgField.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Spray a
zonalPol n lambda = 
  jackPol n lambda (fromInteger 2) 'C'

-- | Symbolic Schur polynomial
schurPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> QSpray 
schurPol' = schurPol

-- | Symbolic Schur polynomial
schurPol :: forall a. (Eq a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Spray a
schurPol n lambda 
  | n < 0 = error "schurPol: negative number of variables."
  | not (_isPartition lambda) = 
      error "schurPol: invalid integer partition."
  | n == 0 = if null lambda then unitSpray else zeroSpray
  | otherwise = sch n 1 lambda arr0
      where
        x = map lone [1 .. n] :: [Spray a]
        nll = _N lambda lambda
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: 
          Int -> Int -> [Int] -> Array (Int,Int) (Maybe (Spray a)) -> Spray a
        sch m k nu arr
          | null nu || nu!!0 == 0 || m == 0 = unitSpray
          | length nu > m && nu !! m > 0 = zeroSpray
          | m == 1 = x !! 0 ^**^ nu !! 0
          | isJust (arr ! (_N lambda nu, m)) = 
              fromJust $ arr ! (_N lambda nu, m)
          | otherwise = s
            where
              s = go (sch (m-1) 1 nu arr) k
              go :: Spray a -> Int -> Spray a
              go !ss ii
                | length nu < ii || nu!!(ii-1) == 0 = ss
                | otherwise =
                  let u = nu!!(ii-1) in
                  if length nu == ii && u > 0 || u > nu !! ii
                    then
                      let nu' = (element (ii-1) .~ u-1) nu in
                      if u > 1
                        then
                          go (ss ^+^ ((x!!(m-1)) ^*^ sch m ii nu' arr)) 
                              (ii + 1)
                        else
                          if nu' !! 0 == 0
                            then
                              go (ss ^+^ (x!!(m-1))) (ii + 1)
                            else
                              let arr' = 
                                    arr // [((_N lambda nu, m), Just ss)] in
                              go (ss ^+^ ((x!!(m-1)) ^*^ sch (m-1) 1 nu' arr')) 
                                  (ii + 1)
                    else
                      go ss (ii + 1)

-- | Symbolic skew Schur polynomial
skewSchurPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> QSpray
skewSchurPol' = skewSchurPol

-- | Symbolic skew Schur polynomial
skewSchurPol :: forall a. (Eq a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Spray a
skewSchurPol n lambda mu =
  case isSkewPartition lambda mu of
    False -> error "skewSchurPol: invalid skew partition."
    True  -> DM.foldlWithKey' f zeroSpray lrCoefficients
  where
    lrCoefficients = skewSchurLRCoefficients lambda mu
    f :: Spray a -> Partition -> Int -> Spray a
    f spray nu k = spray ^+^ (_fromInt k) *^ (schurPol n nu)
