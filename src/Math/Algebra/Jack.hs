{-|
Module      : Math.Algebra.JackPol
Description : Evaluation of Jack polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

Evaluation of Jack polynomials, zonal polynomials, and Schur polynomials. 
See README for examples and references.
-}

{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.Jack
  (jack, zonal, schur, skewSchur)
  where
import qualified Algebra.Additive as AA
import qualified Algebra.Ring     as AR
import Control.Lens               ( (.~), element )
import Data.Array                 ( Array, (!), (//), listArray )
import Data.Maybe                 ( fromJust, isJust )
import qualified Data.Map.Strict  as DM
import Math.Algebra.Jack.Internal ( _N, hookLengths
                                  , _betaratio, _isPartition
                                  , Partition, skewSchurLRCoefficients )
import Numeric.SpecFunctions      ( factorial )

-- | Evaluation of Jack polynomial
jack :: forall a. (Fractional a, Ord a) 
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers
  -> a         -- ^ alpha parameter
  -> a
jack []       _      _     = error "jack: empty list of variables"
jack x@(x0:_) lambda alpha =
  case _isPartition lambda && alpha > 0 of
    False -> if _isPartition lambda
      then error "jack: alpha must be strictly positive"
      else error "jack: invalid integer partition"
    True -> jac (length x) 0 lambda lambda arr0 1
      where
      nll = _N lambda lambda
      n = length x
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then 1
        else product $ map (\i -> alpha * fromIntegral i + 1) [1 .. nu0-1]
      jac :: Int -> Int -> [Int] -> [Int] -> Array (Int,Int) (Maybe a) -> a -> a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = 1
        | length nu > m && nu!!m > 0 = 0
        | m == 1 = x0 ^ (nu!!0) * theproduct (nu!!0)
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (jac (m-1) 0 nu nu arr 1 * beta * x!!(m-1) ^ (sum mu - sum nu))
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
                            go (ss + gamma * x!!(m-1)^ sum mu) (ii + 1)
                          else
                            let arr' = arr // [((_N lambda nu, m), Just ss)] in
                            let jck = jac (m-1) 0 nu' nu' arr' 1 in
                            let jck' = jck * gamma *
                                        x!!(m-1) ^ (sum mu - sum nu') in
                            go (ss+jck') (ii+1)
                  else
                    go ss (ii+1)

-- | Evaluation of zonal polynomial
zonal :: (Fractional a, Ord a) 
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers
  -> a
zonal x lambda = c * jck
  where
    k = sum lambda
    jlambda = product (hookLengths lambda 2)
    c = 2^k * realToFrac (factorial k) / jlambda
    jck = jack x lambda 2

-- | Evaluation of Schur polynomial
schur :: forall a. AR.C a 
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
          | null nu || nu!!0 == 0 || m == 0 = AR.one
          | length nu > m && nu!!m > 0 = AA.zero
          | m == 1 = AR.product (replicate (nu!!0) x0)
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
                          go (ss AA.+ x!!(m-1) AR.* sch m ii nu' arr) (ii + 1)
                        else
                          if nu' !! 0 == 0
                            then
                              go (ss AA.+ x!!(m-1)) (ii + 1)
                            else
                              let arr' = arr // [((_N lambda nu, m), Just ss)] in
                              go (ss AA.+ x!!(m-1) AR.* sch (m-1) 1 nu' arr') (ii + 1)
                    else
                      go ss (ii+1)

-- | Evaluation of a skew Schur polynomial
skewSchur :: forall a. AR.C a 
  => [a]       -- ^ values of the variables
  -> Partition -- ^ partition of integers, the outer partition 
  -> Partition -- ^ partition of integers, the inner partition
  -> a
skewSchur xs lambda mu = DM.foldlWithKey' f AA.zero lrCoefficients
  where
    lrCoefficients = skewSchurLRCoefficients lambda mu
    f :: a -> Partition -> a -> a
    f x nu k = x AA.+ k AR.* (schur xs nu)

