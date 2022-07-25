{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Jack4
  (schur, jack, zonal)
  where
import Control.Lens                             ( (.~), element )
import Data.Array                               ( Array, (!), (//), listArray )
import Data.Maybe                               ( fromJust, isJust )
import Internal                                 ( _N, hookLengths, _betaratio, _isPartition )
import Numeric.SpecFunctions                    ( factorial )

jack :: forall a. (Fractional a, Ord a) => [a] -> [Int] -> a -> a
jack x lambda alpha =
  case _isPartition lambda && alpha > 0 of
    False -> if _isPartition lambda
      then error "alpha must be strictly positive"
      else error "lambda is not a valid integer partition"
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
        | null nu || head nu == 0 || m == 0 = 1
        | length nu > m && nu!!m > 0 = 0
        | m == 1 = head x ^ head nu * theproduct (head nu)
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
                        if head nu' == 0
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

zonal :: (Fractional a, Ord a) => [a] -> [Int] -> a
zonal x lambda = c * jck
  where
    k = sum lambda
    jlambda = product (hookLengths lambda 2)
    c = 2^k * realToFrac (factorial k) / jlambda
    jck = jack x lambda 2

schur :: forall a. Fractional a => [a] -> [Int] -> a
schur x lambda =
  case _isPartition lambda of
    False -> error "lambda is not a valid integer partition"
    True -> sch n 1 lambda arr0
      where
        nll = _N lambda lambda
        n = length x
        arr0 = listArray ((1, 1), (nll, n)) (replicate (nll * n) Nothing)
        sch :: Int -> Int -> [Int] -> Array (Int,Int) (Maybe a) -> a
        sch m k nu arr
          | null nu || head nu == 0 || m == 0 = 1
          | length nu > m && nu!!m > 0 = 0
          | m == 1 = head x ^ head nu
          | isJust (arr ! (_N lambda nu, m)) = fromJust $ arr ! (_N lambda nu, m)
          | otherwise = s
            where
              s = go (sch (m-1) 1 nu arr) k
              go :: Fractional a => a -> Int -> a
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
                          if head nu' == 0
                            then
                              go (ss + x!!(m-1)) (ii + 1)
                            else
                              let arr' = arr // [((_N lambda nu, m), Just ss)] in
                              go (ss + x!!(m-1) * sch (m-1) 1 nu' arr') (ii + 1)
                    else
                      go ss (ii+1)
