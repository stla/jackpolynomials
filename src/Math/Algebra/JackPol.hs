{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.JackPol
  (schurPol, jackPol, zonalPol)
  where
import qualified Algebra.Ring               as AR
import           Control.Lens               ( (.~), element )
import           Data.Array                 ( Array, (!), (//), listArray )
import           Data.Maybe                 ( fromJust, isJust )
import           Math.Algebra.Jack.Internal ( _betaratio, hookLengths, _N
                                            , _isPartition, Partition )
import           Math.Algebra.Hspray        ( (*^), (^**^), (^*^), (^+^)
                                            , lone, Spray
                                            , zeroSpray, unitSpray )
import           Numeric.SpecFunctions      ( factorial )

-- | Symbolic Jack polynomial
jackPol :: forall a. (Fractional a, Ord a, AR.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> a         -- ^ alpha parameter
  -> Spray a
jackPol n lambda alpha =
  case _isPartition lambda && alpha > 0 of
    False -> if _isPartition lambda
      then error "jackPol: alpha must be strictly positive"
      else error "jackPol: invalid integer partition"
    True -> jac (length x) 0 lambda lambda arr0 1
      where
      nll = _N lambda lambda
      x = map lone [1 .. n] :: [Spray a]
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then 1
        else product $ map (\i -> alpha * fromIntegral i + 1) [1 .. nu0-1]
      jac :: Int -> Int -> Partition -> Partition -> Array (Int,Int) (Maybe (Spray a)) -> a -> Spray a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = unitSpray
        | length nu > m && nu!!m > 0 = zeroSpray
        | m == 1 = theproduct (nu!!0) *^ (x!!0 ^**^ nu!!0) 
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (beta *^ (jac (m-1) 0 nu nu arr 1 ^*^ ((x!!(m-1)) ^**^ (sum mu - sum nu))))
                (max 1 k)
            go :: Spray a -> Int -> Spray a
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
                        go (ss ^+^ jac m ii mu nu' arr gamma) (ii + 1)
                      else
                        if nu'!!0 == 0
                          then
                            go (ss ^+^ (gamma *^ (x!!(m-1) ^**^ sum mu))) (ii + 1)
                          else
                            let arr' = arr // [((_N lambda nu, m), Just ss)] in
                            let jck = jac (m-1) 0 nu' nu' arr' 1 in
                            let jck' = gamma *^ (jck ^*^ 
                                        (x!!(m-1) ^**^ (sum mu - sum nu'))) in
                            go (ss ^+^ jck') (ii+1)
                  else
                    go ss (ii+1)

-- | Symbolic zonal polynomial
zonalPol :: (Fractional a, Ord a, AR.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Spray a
zonalPol n lambda = c *^ jck
  where
    k = sum lambda
    jlambda = product (hookLengths lambda 2)
    c = 2^k * realToFrac (factorial k) / jlambda
    jck = jackPol n lambda 2

-- | Symbolic Schur polynomial
schurPol :: forall a. (Real a, Ord a, AR.C a)
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
