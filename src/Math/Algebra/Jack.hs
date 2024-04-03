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
  (jack', zonal', schur', skewSchur', jack, zonal, schur, skewSchur)
  where
import           Prelude 
  hiding ((*), (+), (-), (/), (^), (*>), product, sum, fromIntegral, fromInteger)
import           Algebra.Additive           ( (+), (-), sum, zero )
import           Algebra.Ring               ( (*), product, one, (^), fromInteger )
import           Algebra.ToInteger          ( fromIntegral ) 
import qualified Algebra.Field              as AlgField
import qualified Algebra.Ring               as AlgRing
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import           Data.Maybe                 ( fromJust, isJust )
import qualified Data.Map.Strict            as DM
import           Math.Algebra.Jack.Internal ( (.^), _N, jackCoeffC
                                            , jackCoeffP, jackCoeffQ
                                            , _betaratio, _isPartition
                                            , Partition, skewSchurLRCoefficients
                                            , isSkewPartition, _fromInt )

-- | Evaluation of Jack polynomial
jack' 
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ partition of integers
  -> Rational   -- ^ Jack parameter
  -> Char       -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> Rational
jack' = jack

-- | Evaluation of Jack polynomial
jack :: forall a. AlgField.C a
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers
  -> a         -- ^ Jack parameter
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> a
jack []       _      _     _     = error "jack: empty list of variables"
jack x@(x0:_) lambda alpha which =
  case _isPartition lambda of
    False -> error "jack: invalid integer partition"
    True -> case which of 
      'J' -> resultJ
      'C' -> jackCoeffC lambda alpha * resultJ
      'P' -> jackCoeffP lambda alpha * resultJ
      'Q' -> jackCoeffQ lambda alpha * resultJ
      _   -> error "jack: please use 'J', 'C', 'P' or 'Q' for last argument"
      where
      resultJ = jac (length x) 0 lambda lambda arr0 one
      nll = _N lambda lambda
      n = length x
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then one
        else product $ map (\i -> one + i .^ alpha) [1 .. nu0-1]
      jac :: Int -> Int -> [Int] -> [Int] -> Array (Int,Int) (Maybe a) -> a -> a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = one
        | length nu > m && nu!!m > 0      = zero
        | m == 1 = x0 ^ (fromIntegral $ nu!!0) * theproduct (nu!!0)
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (jac (m-1) 0 nu nu arr one * beta * x!!(m-1) ^ (fromIntegral $ sum mu - sum nu))
                (max 1 k)
            go :: a -> Int -> a
            go !ss ii
              | length nu < ii || nu!!(ii-1) == 0 = ss
              | otherwise =
                let u = nu!!(ii-1) in
                if length nu == ii && u > 0 || u > nu!!ii
                  then
                    let nu' = (element (ii-1) .~ u-1) nu in
                    let gamma = beta * _betaratio mu nu ii alpha in
                    if u > 1
                      then
                        go (ss + jac m ii mu nu' arr gamma) (ii + 1)
                      else
                        if nu' !! 0 == 0
                          then
                            go (ss + gamma * x!!(m-1)^ (fromIntegral $ sum mu)) (ii + 1)
                          else
                            let arr' = arr // [((_N lambda nu, m), Just ss)] in
                            let jck  = jac (m-1) 0 nu' nu' arr' one in
                            let jck' = jck * gamma *
                                        x!!(m-1) ^ (fromIntegral $ sum mu - sum nu') in
                            go (ss + jck') (ii + 1)
                  else
                    go ss (ii + 1)

-- | Evaluation of zonal polynomial
zonal' 
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ partition of integers
  -> Rational
zonal' = zonal

-- | Evaluation of zonal polynomial
zonal :: AlgField.C a
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
        nll = _N lambda lambda
        n = length x
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: Int -> Int -> [Int] -> Array (Int,Int) (Maybe a) -> a
        sch m k nu arr
          | null nu || nu!!0 == 0 || m == 0 = one
          | length nu > m && nu!!m > 0 = zero
          | m == 1 = product (replicate (nu!!0) x0)
          | isJust (arr ! (_N lambda nu, m)) = fromJust $ arr ! (_N lambda nu, m)
          | otherwise = s
            where
              s = go (sch (m-1) 1 nu arr) k
              go :: a -> Int -> a
              go !ss ii
                | length nu < ii || nu!!(ii-1) == 0 = ss
                | otherwise =
                  let u = nu!!(ii-1) in
                  if length nu == ii && u > 0 || u > nu !! ii
                    then
                      let nu' = (element (ii-1) .~ u-1) nu in
                      if u > 1
                        then
                          go (ss + x!!(m-1) * sch m ii nu' arr) (ii + 1)
                        else
                          if nu' !! 0 == 0
                            then
                              go (ss + x!!(m-1)) (ii + 1)
                            else
                              let arr' = arr // [((_N lambda nu, m), Just ss)] in
                              go (ss + x!!(m-1) * sch (m-1) 1 nu' arr') (ii + 1)
                    else
                      go ss (ii + 1)

-- | Evaluation of a skew Schur polynomial
skewSchur' 
  :: [Rational] -- ^ values of the variables
  -> Partition  -- ^ the outer partition of the skew partition
  -> Partition  -- ^ the inner partition of the skew partition
  -> Rational
skewSchur' = skewSchur

-- | Evaluation of a skew Schur polynomial
skewSchur :: forall a. AlgRing.C a 
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

