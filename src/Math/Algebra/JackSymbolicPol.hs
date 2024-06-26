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
  ( jackSymbolicPol, jackSymbolicPol', skewJackSymbolicPol, skewJackSymbolicPol' )
  where
import           Prelude 
  hiding ((/), (^), (*>), product, fromIntegral, fromInteger, recip)
import           Algebra.Field              ( recip )
import qualified Algebra.Field              as AlgField
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import qualified Data.HashMap.Strict        as HM
import qualified Data.Map.Strict            as DM
import           Data.Maybe                 ( fromJust, isJust )
import           Math.Algebra.Jack.Internal ( _betaRatioOfSprays
                                            , jackSymbolicCoeffC
                                            , jackSymbolicCoeffPinv
                                            , jackSymbolicCoeffQinv
                                            , _N, _isPartition, Partition
                                            , isSkewPartition
                                            , skewSymbolicJackInMSPbasis )
import           Math.Algebra.Hspray        ( FunctionLike (..), (.^)
                                            , lone, lone'
                                            , ParametricSpray, ParametricQSpray
                                            , Spray, asRatioOfSprays
                                            , zeroSpray, unitSpray
                                            , productOfSprays
                                            , fromList )
import           Math.Combinat.Permutations ( permuteMultiset )


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
    False -> error "jackSymbolicPol: invalid integer partition."
    True -> case which of 
      'J' -> resultJ
      'C' -> jackSymbolicCoeffC lambda *^ resultJ
      'P' -> recip (asRatioOfSprays (jackSymbolicCoeffPinv lambda)) *^ resultJ 
      'Q' -> recip (asRatioOfSprays (jackSymbolicCoeffQinv lambda)) *^ resultJ
      _   -> error 
        "jackSymbolicPol: please use 'J', 'C', 'P' or 'Q' for last argument."
      where
      alpha = lone 1 :: Spray a
      jck m kappa arr = jac m 0 kappa kappa arr 
      resultJ = jck n lambda arr0 
      nll = _N lambda lambda
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      jac :: Int -> Int -> Partition -> Partition 
             -> Array (Int,Int) (Maybe (ParametricSpray a)) 
             -> ParametricSpray a
      jac m k mu nu arr 
        | null nu || nu0 == 0 || m == 0 = unitSpray
        | ellNu > m && nu !! m > 0      = zeroSpray
        | m == 1                        = 
            if nu0 == 1
              then 
                lone 1
              else 
                let sprays = [i .^ alpha ^+^ unitSpray | i <- [1 .. nu0-1]] in
                asRatioOfSprays (productOfSprays sprays) *^ x nu0
        | k == 0 && isJust maybe_pspray =
            fromJust $ maybe_pspray
        | otherwise = s
          where
            nu0 = nu !! 0
            ellNu = length nu
            x = lone' m
            _N_lambda_nu_m = (_N lambda nu, m)
            maybe_pspray = arr ! _N_lambda_nu_m
            wMu = sum mu
            jck' kappa array = jck (m-1) kappa array ^*^ x (wMu - sum kappa)
            s = go (jck' nu arr) (max 1 k)
            go :: ParametricSpray a -> Int -> ParametricSpray a
            go !ss ii
              | ellNu < ii || u == 0 = 
                  ss
              | ellNu == ii && u > 0 || u > nu !! ii = 
                  go (ss ^+^ tt) (ii + 1)
              | otherwise = 
                  go ss (ii + 1)
                where
                  jj = ii - 1
                  u = nu !! jj
                  nu' = (element jj .~ u - 1) nu
                  gamma = _betaRatioOfSprays mu nu ii
                  tt = gamma *^ pspray
                    where
                      pspray
                        | u > 1 =
                            jac m ii mu nu' arr 
                        | nu' !! 0 == 0 =
                            x wMu
                        | otherwise =
                            jck' nu' (arr // [(_N_lambda_nu_m, Just ss)]) 

skewJackSymbolicPol :: 
    (Eq a, AlgField.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Char      -- ^ which skew Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> ParametricSpray a
skewJackSymbolicPol n lambda mu which 
  | n < 0 = 
      error "skewJackSymbolicPol: negative number of variables."
  | not (isSkewPartition lambda mu) = 
      error "skewJackSymbolicPol: invalid skew partition."
  | not (which `elem` ['J', 'C', 'P', 'Q']) = 
      error "skewJackSymbolicPol: please use 'J', 'C', 'P' or 'Q' for last argument."
  | n == 0 = 
      if lambda == mu then unitSpray else zeroSpray
  | otherwise =
      HM.unions sprays
  where
    msCombo = 
      DM.filter
        ((<= n) . fst)
          (DM.mapWithKey 
            (\kappa rOS -> (length kappa, rOS)) 
              (skewSymbolicJackInMSPbasis which lambda mu))
    sprays = 
      map (
        \(kappa, (l, rOS)) -> 
          fromList
            (zip 
              (permuteMultiset (kappa ++ replicate (n - l) 0)) 
              (repeat rOS))
        ) (DM.assocs msCombo)

skewJackSymbolicPol' :: 
     Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Char      -- ^ which skew Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> ParametricQSpray 
skewJackSymbolicPol' = skewJackSymbolicPol