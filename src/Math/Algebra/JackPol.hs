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
  (jackPol, zonalPol, schurPol, skewSchurPol)
  where
import Prelude hiding ((*), (+), (-), (/), (^), (*>), product, sum, fromIntegral, fromInteger)
import           Algebra.Additive           
import           Algebra.Module             
import           Algebra.Ring
import           Algebra.ToInteger           
import qualified Algebra.Module             as AlgMod
import qualified Algebra.Field              as AlgField
import qualified Algebra.Ring               as AlgRing
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import qualified Data.Map.Strict            as DM
import           Data.Maybe                 ( fromJust, isJust )
import           Math.Algebra.Jack.Internal ( (.^), _betaratio, _productHookLengths
                                            , _N, _isPartition, Partition
                                            , jackCoeffP, jackCoeffQ
                                            , skewSchurLRCoefficients
                                            , isSkewPartition, _fromInt )
import           Math.Algebra.Hspray        ( (*^), (^**^), (^*^), (^+^)
                                            , lone, Spray
                                            , zeroSpray, unitSpray )

-- | Symbolic Jack polynomial
jackPol :: forall a. (AlgField.C a, Ord a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> a         -- ^ alpha parameter
  -> String    -- ^ which Jack polynomial, @"J"@, @"P"@ or @"Q"@
  -> Spray a
jackPol n lambda alpha which =
  case _isPartition lambda && alpha > zero of
    False -> if _isPartition lambda
      then error "jackPol: alpha must be strictly positive"
      else error "jackPol: invalid integer partition"
    True -> case which of 
      "J" -> resultJ
      "P" -> jackCoeffP lambda alpha *> resultJ
      "Q" -> jackCoeffQ lambda alpha *> resultJ
      _   -> error "jackPol: please use \"J\", \"P\" or \"Q\" for last argument"
      where
      resultJ = jac (length x) 0 lambda lambda arr0 one
      nll = _N lambda lambda
      x = map lone [1 .. n] :: [Spray a]
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then one
        else product $ map (\i -> i .^ alpha + one) [1 .. nu0-1]
      jac :: Int -> Int -> Partition -> Partition -> Array (Int,Int) (Maybe (Spray a)) -> a -> Spray a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = unitSpray
        | length nu > m && nu!!m > 0      = zeroSpray
        | m == 1                          = theproduct (nu!!0) *^ (x!!0 ^**^ nu!!0) 
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (beta *^ (jac (m-1) 0 nu nu arr one ^*^ ((x!!(m-1)) ^**^ (sum mu - sum nu))))
                (max 1 k)
            go :: Spray a -> Int -> Spray a
            go !ss ii
              | length nu < ii || nu!!(ii-1) == 0 = ss
              | otherwise =
                let u = nu!!(ii-1) in
                if length nu == ii && u > 0 || u > nu!!ii
                  then
                    let nu'   = (element (ii-1) .~ u-1) nu in
                    let gamma = beta * _betaratio mu nu ii alpha in
                    if u > 1
                      then
                        go (ss ^+^ jac m ii mu nu' arr gamma) (ii + 1)
                      else
                        if nu'!!0 == 0
                          then
                            go (ss ^+^ (gamma *^ (x!!(m-1) ^**^ sum mu))) (ii + 1)
                          else
                            let arr' = arr // [((_N lambda nu, m), Just ss)] in
                            let jck  = jac (m-1) 0 nu' nu' arr' one in
                            let jck' = gamma *^ (jck ^*^ 
                                        (x!!(m-1) ^**^ (sum mu - sum nu'))) in
                            go (ss ^+^ jck') (ii + 1)
                  else
                    go ss (ii + 1)

-- | Symbolic zonal polynomial
zonalPol :: forall a. (AlgField.C a, Ord a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Spray a
zonalPol n lambda = fromInteger c * (AlgField.recip jlambda *> jck)
  where
    k = fromIntegral $ sum lambda
    jlambda = _productHookLengths lambda (fromInteger 2) :: a
    c = 2^k * (product [2 .. k])
    jck = jackPol n lambda (fromInteger 2) "J"

-- | Symbolic Schur polynomial
schurPol :: forall a. (Ord a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Spray a
schurPol n lambda =
  case _isPartition lambda of
    False -> error "schurPol: invalid integer partition"
    True -> sch n 1 lambda arr0
      where
        x = map lone [1 .. n] :: [Spray a]
        nll = _N lambda lambda
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: Int -> Int -> [Int] -> Array (Int,Int) (Maybe (Spray a)) -> Spray a
        sch m k nu arr
          | null nu || nu!!0 == 0 || m == 0 = unitSpray
          | length nu > m && nu!!m > 0 = zeroSpray
          | m == 1 = x!!0 ^**^ nu!!0
          | isJust (arr ! (_N lambda nu, m)) = fromJust $ arr ! (_N lambda nu, m)
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
                          go (ss ^+^ ((x!!(m-1)) ^*^ sch m ii nu' arr)) (ii + 1)
                        else
                          if nu'!!0 == 0
                            then
                              go (ss ^+^ (x!!(m-1))) (ii + 1)
                            else
                              let arr' = arr // [((_N lambda nu, m), Just ss)] in
                              go (ss ^+^ ((x!!(m-1)) ^*^ sch (m-1) 1 nu' arr')) (ii + 1)
                    else
                      go ss (ii+1)

-- | Symbolic skew Schur polynomial
skewSchurPol :: forall a. (Ord a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Spray a
skewSchurPol n lambda mu =
  case isSkewPartition lambda mu of
    False -> error "skewSchurPol: invalid skew partition"
    True  -> DM.foldlWithKey' f zeroSpray lrCoefficients
  where
    lrCoefficients = skewSchurLRCoefficients lambda mu
    f :: Spray a -> Partition -> Int -> Spray a
    f spray nu k = spray ^+^ (_fromInt' k) AlgMod.*> (schurPol n nu)
    _fromInt' :: Int -> a
    _fromInt' = _fromInt

