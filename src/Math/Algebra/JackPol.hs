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
                                            , lone, lone', Spray, QSpray
                                            , zeroSpray, unitSpray
                                            , fromList )
import           Math.Combinat.Permutations ( permuteMultiset )

-- | Jack polynomial
jackPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Rational  -- ^ Jack parameter
  -> Char      -- ^ which Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> QSpray
jackPol' = jackPol

-- | Jack polynomial
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
      jck m kappa arr = jac m 0 kappa kappa arr 
      resultJ = jck n lambda arr0 
      nll = _N lambda lambda
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      jac :: Int -> Int -> Partition -> Partition 
             -> Array (Int,Int) (Maybe (Spray a)) 
             -> Spray a
      jac m k mu nu arr 
        | null nu || nu0 == 0 || m == 0 = unitSpray
        | ellNu > m && nu !! m > 0      = zeroSpray
        | m == 1                        = 
            if nu0 == 1
              then 
                lone 1
              else 
                let as = [i .^ alpha + one | i <- [1 .. nu0-1]] in
                product as *^ x nu0
        | k == 0 && isJust maybe_spray =
            fromJust $ maybe_spray
        | otherwise = s
          where
            nu0 = nu !! 0
            ellNu = length nu
            x = lone' m
            _N_lambda_nu_m = (_N lambda nu, m)
            maybe_spray = arr ! _N_lambda_nu_m
            wMu = sum mu
            jck' kappa array = jck (m-1) kappa array ^*^ x (wMu - sum kappa)
            s = go (jck' nu arr) (max 1 k)
            go :: Spray a -> Int -> Spray a
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
                  gamma = _betaratio mu nu ii alpha
                  tt = gamma *^ spray
                    where
                      spray
                        | u > 1 =
                            jac m ii mu nu' arr 
                        | nu' !! 0 == 0 =
                            x wMu
                        | otherwise =
                            jck' nu' (arr // [(_N_lambda_nu_m, Just ss)]) 

skewJackPol :: 
    (Eq a, AlgField.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> a         -- ^ Jack parameter
  -> Char      -- ^ which skew Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> Spray a
skewJackPol n lambda mu alpha which 
  | n < 0 = 
      error "skewJackPol: negative number of variables."
  | not (isSkewPartition lambda mu) = 
      error "skewJackPol: invalid skew partition."
  | not (which `elem` ['J', 'C', 'P', 'Q']) = 
      error "skewJackPol: please use 'J', 'C', 'P' or 'Q' for last argument."
  | n == 0 = 
      if lambda == mu then unitSpray else zeroSpray
  | otherwise =
      HM.unions sprays
  where
    msCombo = 
      DM.filter
        ((<= n) . fst)
          (DM.mapWithKey 
            (\kappa coeff -> (length kappa, coeff)) 
              (skewJackInMSPbasis alpha which lambda mu))
    sprays = 
      map (
        \(kappa, (l, coeff)) -> 
          fromList
            (zip 
              (permuteMultiset (kappa ++ replicate (n - l) 0)) 
              (repeat coeff))
        ) (DM.assocs msCombo)

skewJackPol' :: 
     Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Rational  -- ^ Jack parameter
  -> Char      -- ^ which skew Jack polynomial, @'J'@, @'C'@, @'P'@ or @'Q'@
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
        nll = _N lambda lambda
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: 
          Int -> Int -> [Int] -> Array (Int,Int) (Maybe (Spray a)) -> Spray a
        sch m k nu arr
          | null nu || nu0 == 0 || m == 0 = unitSpray
          | ellNu > m && nu !! m > 0 = zeroSpray
          | m == 1 = lone' 1 nu0
          | isJust maybe_spray = 
              fromJust maybe_spray
          | otherwise = s
            where
              nu0 = nu !! 0
              ellNu = length nu
              x = lone m
              _N_lambda_nu_m = (_N lambda nu, m)
              maybe_spray = arr ! _N_lambda_nu_m
              sch' kappa array = sch (m-1) 1 kappa array 
              s = go (sch' nu arr) k
              go :: Spray a -> Int -> Spray a
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
                    tt 
                      | u > 1 =
                          x ^*^ sch m ii nu' arr
                      | nu' !! 0 == 0 =
                          x 
                      | otherwise =
                          x ^*^ sch' nu' (arr // [(_N_lambda_nu_m, Just ss)]) 

-- | Skew Schur polynomial
skewSchurPol' 
  :: Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> QSpray
skewSchurPol' = skewSchurPol

-- | Skew Schur polynomial
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
