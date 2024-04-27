{-|
Module      : Math.Algebra.JackSymbolicPol
Description : Jack polynomials with symbolic Jack parameter.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

Computation of Jack polynomials with a symbolic Jack parameter. 
See README for examples and references.
-}

{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.JackSymbolicPol
  (jackSymbolicPol', jackSymbolicPol)
  where
import           Prelude 
  hiding ((*), (+), (-), (/), (^), (*>), product, sum, fromIntegral, fromInteger, recip)
import           Algebra.Additive           ( (+), (-), sum )
import           Algebra.Ring               ( (*), product, one )
import           Algebra.Field              ( recip )
import qualified Algebra.Field              as AlgField
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import           Data.Maybe                 ( fromJust, isJust )
import           Math.Algebra.Jack.Internal ( _betaRatioOfSprays
                                            , jackSymbolicCoeffC
                                            , jackSymbolicCoeffPinv
                                            , jackSymbolicCoeffQinv
                                            , _N, _isPartition, Partition )
import           Math.Algebra.Hspray        ( (*^), (^**^), (^*^), (^+^), (.^)
                                            , lone, ParametricSpray, ParametricQSpray
                                            , Spray, asRatioOfSprays
                                            , RatioOfSprays, unitRatioOfSprays
                                            , zeroSpray, unitSpray )

-- | Jack polynomial with symbolic Jack parameter
jackSymbolicPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> ParametricQSpray
jackSymbolicPol' = jackSymbolicPol

-- | Jack polynomial with symbolic Jack parameter
jackSymbolicPol :: forall a. (Eq a, AlgField.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> ParametricSpray a
jackSymbolicPol n lambda which =
  case _isPartition lambda of
    False -> error "jackSymbolicPol: invalid integer partition"
    True -> case which of 
      'J' -> resultJ
      'C' -> jackSymbolicCoeffC lambda *^ resultJ
      'P' -> recip (asRatioOfSprays (jackSymbolicCoeffPinv lambda)) *^ resultJ 
      'Q' -> recip (asRatioOfSprays (jackSymbolicCoeffQinv lambda)) *^ resultJ
      _   -> error 
        "jackSymbolicPol: please use 'J', 'C', 'P' or 'Q' for last argument"
      where
      alpha = lone 1 :: Spray a
      resultJ = jac (length x) 0 lambda lambda arr0 one
      nll = _N lambda lambda
      x = map lone [1 .. n] :: [ParametricSpray a]
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> RatioOfSprays a
      theproduct nu0 = if nu0 <= 1
        then unitRatioOfSprays
        else asRatioOfSprays $ product $ map 
              (\i -> i .^ alpha ^+^ unitSpray) 
              [1 .. nu0-1]
      jac :: Int -> Int -> Partition -> Partition 
             -> Array (Int,Int) (Maybe (ParametricSpray a)) 
             -> RatioOfSprays a -> ParametricSpray a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = unitSpray
        | length nu > m && nu!!m > 0      = zeroSpray
        | m == 1                          = theproduct (nu!!0) *^ (x!!0 ^**^ nu!!0) 
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (beta *^ (jac (m-1) 0 nu nu arr unitRatioOfSprays ^*^ ((x!!(m-1)) ^**^ (sum mu - sum nu))))
                (max 1 k)
            go :: ParametricSpray a -> Int -> ParametricSpray a
            go !ss ii
              | length nu < ii || nu!!(ii-1) == 0 = ss
              | otherwise =
                let u = nu!!(ii-1) in
                if length nu == ii && u > 0 || u > nu!!ii
                  then
                    let nu'   = (element (ii-1) .~ u-1) nu in
                    let gamma = _betaRatioOfSprays mu nu ii * beta in
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

