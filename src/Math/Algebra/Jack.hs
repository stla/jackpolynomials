{-|
Module      : Math.Algebra.Jack
Description : Evaluation of Jack polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

Evaluation of Jack polynomials, zonal polynomials, Schur polynomials and skew Schur polynomials. 
See README for examples and references.
-}

{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.Jack
  (Partition, jack', zonal', schur', skewSchur', jack, zonal, schur, skewSchur)
  where
import           Prelude 
  hiding ((*), (+), (-), (/), (^), (*>), product, fromIntegral, fromInteger)
import           Algebra.Additive           ( (+), (-), zero )
import           Algebra.Ring               ( (*), product, one, (^), fromInteger )
import qualified Algebra.Field              as AlgField
import qualified Algebra.Ring               as AlgRing
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import           Data.Maybe                 ( fromJust, isJust )
import qualified Data.Map.Strict            as DM
import           Math.Algebra.Jack.Internal ( _N, jackCoeffC
                                            , jackCoeffP, jackCoeffQ
                                            , _betaratio, _isPartition
                                            , Partition, skewSchurLRCoefficients
                                            , isSkewPartition, _fromInt )
import Math.Algebra.Hspray                  ( (.^) )

-- | Evaluation of Jack polynomial
jack' 
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ partition of integers
  -> Rational   -- ^ Jack parameter
  -> Char       -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> Rational
jack' = jack

-- | Evaluation of Jack polynomial
jack :: forall a. (Eq a, AlgField.C a)
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers
  -> a         -- ^ Jack parameter
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> a
jack []       _      _     _     = error "jack: empty list of variables."
jack x@(x0:_) lambda alpha which =
  case _isPartition lambda of
    False -> error "jack: invalid integer partition."
    True -> case which of 
      'J' -> resultJ
      'C' -> jackCoeffC lambda alpha * resultJ
      'P' -> jackCoeffP lambda alpha * resultJ
      'Q' -> jackCoeffQ lambda alpha * resultJ
      _   -> error "jack: please use 'J', 'C', 'P' or 'Q' for last argument."
      where
      jck m kappa arr = jac m 0 kappa kappa arr 
      n = length x
      resultJ = jck n lambda arr0 
      nll = _N lambda lambda
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      jac :: Int -> Int -> Partition -> Partition 
             -> Array (Int,Int) (Maybe a) 
             -> a
      jac m k mu nu arr 
        | null nu || nu0 == 0 || m == 0 = one
        | ellNu > m && nu !! m > 0      = zero
        | m == 1                        = 
            if nu0 == 1
              then 
                x0
              else 
                let as = [i .^ alpha + one | i <- [1 .. nu0-1]] in
                product as * x0 ^ (toInteger nu0)
        | k == 0 && isJust maybe_a =
            fromJust $ maybe_a
        | otherwise = s
          where
            nu0 = nu !! 0
            ellNu = length nu
            xm = x !! (m - 1)
            xmi i = xm ^ (toInteger i)
            _N_lambda_nu_m = (_N lambda nu, m)
            maybe_a = arr ! _N_lambda_nu_m
            wMu = sum mu
            jck' kappa array = jck (m-1) kappa array * (xmi (wMu - sum kappa))
            s = go (jck' nu arr) (max 1 k)
            go :: a -> Int -> a
            go !ss ii
              | ellNu < ii || u == 0 = 
                  ss
              | ellNu == ii && u > 0 || u > nu !! ii = 
                  go (ss + tt) (ii + 1)
              | otherwise = 
                  go ss (ii + 1)
                where
                  jj = ii - 1
                  u = nu !! jj
                  nu' = (element jj .~ u - 1) nu
                  gamma = _betaratio mu nu ii alpha
                  tt = gamma * y
                    where
                      y
                        | u > 1 =
                            jac m ii mu nu' arr 
                        | nu' !! 0 == 0 =
                            xmi wMu
                        | otherwise =
                            jck' nu' (arr // [(_N_lambda_nu_m, Just ss)]) 

-- | Evaluation of zonal polynomial
zonal' 
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ partition of integers
  -> Rational
zonal' = zonal

-- | Evaluation of zonal polynomial
zonal :: (Eq a, AlgField.C a)
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers
  -> a
zonal x lambda = jack x lambda (fromInteger 2) 'C'

-- | Evaluation of Schur polynomial
schur'
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ partition of integers 
  -> Rational
schur' = schur

-- | Evaluation of Schur polynomial
schur :: forall a. AlgRing.C a 
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers 
  -> a
schur []       _      = error "schur: empty list of variables"
schur x@(x0:_) lambda =
  case _isPartition lambda of
    False -> error "schur: invalid integer partition"
    True -> sch n 1 lambda arr0
      where
        n = length x
        nll = _N lambda lambda
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: 
          Int -> Int -> [Int] -> Array (Int,Int) (Maybe a) -> a
        sch m k nu arr
          | null nu || nu0 == 0 || m == 0 = one
          | ellNu > m && nu !! m > 0 = zero
          | m == 1 = x0 ^ (toInteger nu0)
          | isJust maybe_a = 
              fromJust maybe_a
          | otherwise = s
            where
              nu0 = nu !! 0
              ellNu = length nu
              xm = x !! (m - 1)
              _N_lambda_nu_m = (_N lambda nu, m)
              maybe_a = arr ! _N_lambda_nu_m
              sch' kappa array = sch (m-1) 1 kappa array 
              s = go (sch' nu arr) k
              go :: a -> Int -> a
              go !ss ii
                | ellNu < ii || u == 0 = 
                    ss
                | ellNu == ii && u > 0 || u > nu !! ii = 
                    go (ss + tt) (ii + 1)
                | otherwise = 
                    go ss (ii + 1)
                  where
                    jj = ii - 1
                    u = nu !! jj
                    nu' = (element jj .~ u - 1) nu
                    tt 
                      | u > 1 =
                          xm * sch m ii nu' arr
                      | nu' !! 0 == 0 =
                          xm 
                      | otherwise =
                          xm * sch' nu' (arr // [(_N_lambda_nu_m, Just ss)]) 

-- | Evaluation of a skew Schur polynomial
skewSchur' 
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ the outer partition of the skew partition
  -> Partition  -- ^ the inner partition of the skew partition
  -> Rational
skewSchur' = skewSchur

-- | Evaluation of a skew Schur polynomial
skewSchur :: forall a. (Eq a, AlgRing.C a) 
  => [a]       -- ^ values of the variables
  -> Partition -- ^ the outer partition of the skew partition
  -> Partition -- ^ the inner partition of the skew partition
  -> a
skewSchur xs lambda mu = 
  if isSkewPartition lambda mu 
    then DM.foldlWithKey' f zero lrCoefficients
    else error "skewSchur: invalid skew partition"
  where
    lrCoefficients = skewSchurLRCoefficients lambda mu
    f :: a -> Partition -> Int -> a
    f x nu k = x + (_fromInt k) * (schur xs nu)

