{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module JackPol
  (schurPol, jackPol, zonalPol)
  where
import qualified Algebra.Ring as AR
import           Control.Lens ( (.~), element )
import           Data.Array   ( Array, (!), (//), listArray )
import           Data.Maybe   ( fromJust, isJust )
import           Internal     ( _betaratio', hookLengths, _N
                              , _isPartition, Partition )
import           MultiPol     ( (*^), (^**^), (^*^), (^+^)
                              , constant, lone, Polynomial )
import Numeric.SpecFunctions  ( factorial )

-- | Symbolic Jack polynomial
jackPol :: forall a. (Fractional a, Ord a, AR.C a) 
  => Int -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> a -- ^ alpha parameter
  -> Polynomial a
jackPol n lambda alpha =
  case _isPartition lambda && alpha > 0 of
    False -> if _isPartition lambda
      then error "alpha must be strictly positive"
      else error "lambda is not a valid integer partition"
    True -> jac (length x) 0 lambda lambda arr0 1
      where
      nll = _N lambda lambda
      x = map lone [1 .. n] :: [Polynomial a]
      arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then AR.one
        else AR.product $ map (\i -> alpha * fromIntegral i + 1) [1 .. nu0-1]
      jac :: Int -> Int -> Partition -> Partition -> Array (Int,Int) (Maybe (Polynomial a)) -> a -> Polynomial a
      jac m k mu nu arr beta
        | null nu || head nu == 0 || m == 0 = constant 1
        | length nu > m && nu!!m > 0 = constant 0
        | m == 1 = theproduct (head nu) *^ (head x ^**^ head nu) 
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
            s = go (beta *^ (jac (m-1) 0 nu nu arr 1 ^*^ ((x!!(m-1)) ^**^ (sum mu - sum nu))))
                (max 1 k)
            go :: Polynomial a -> Int -> Polynomial a
            go !ss ii
              | length nu < ii || nu!!(ii-1) == 0 = ss
              | otherwise =
                let u = nu!!(ii-1) in
                if length nu == ii && u > 0 || u > nu!!ii
                  then
                    let nu' = (element (ii-1) .~ u-1) nu in
                    let gamma = beta * _betaratio' mu nu ii alpha in
                    if u > 1
                      then
                        go (ss ^+^ jac m ii mu nu' arr gamma) (ii + 1)
                      else
                        if head nu' == 0
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
  => Int -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Polynomial a
zonalPol n lambda = c *^ jck
  where
    k = sum lambda
    jlambda = product (hookLengths lambda 2)
    c = 2^k * realToFrac (factorial k) / jlambda
    jck = jackPol n lambda 2

-- | Symbolic Schur polynomial
schurPol :: 
  Int -- ^ number of variables
  -> Partition -- ^ partition of integers
  -> Polynomial Int
schurPol n lambda =
  case _isPartition lambda of
    False -> error "lambda is not a valid integer partition"
    True -> sch n 1 lambda arr0
      where
        x = map lone [1 .. n] :: [Polynomial Int]
        nll = _N lambda lambda
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: Int -> Int -> [Int] -> Array (Int,Int) (Maybe (Polynomial Int)) -> Polynomial Int
        sch m k nu arr
          | null nu || head nu == 0 || m == 0 = constant 1
          | length nu > m && nu!!m > 0 = constant 0
          | m == 1 = head x ^**^ head nu
          | isJust (arr ! (_N lambda nu, m)) = fromJust $ arr ! (_N lambda nu, m)
          | otherwise = s
            where
              s = go (sch (m-1) 1 nu arr) k
              go :: Polynomial Int -> Int -> Polynomial Int
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
                          if head nu' == 0
                            then
                              go (ss ^+^ (x!!(m-1))) (ii + 1)
                            else
                              let arr' = arr // [((_N lambda nu, m), Just ss)] in
                              go (ss ^+^ ((x!!(m-1)) ^*^ sch (m-1) 1 nu' arr')) (ii + 1)
                    else
                      go ss (ii+1)
