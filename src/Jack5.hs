{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Jack5
  (schur, jack, zonal)
  where
import           Control.Lens                             hiding (iconcatMap)
import           Data.Array.Unboxed
import           Data.List.Index                          (iconcatMap)
import           Data.Maybe
import           Math.Combinat.Partitions.Integer.IntList (_dualPartition,
                                                           _isPartition)
import           Numeric.SpecFunctions                    (factorial)

_ij :: [Int] -> ([Int],[Int])
_ij lambda =
  (
    iconcatMap (\i a ->  replicate a (i+1)) lambda,
    concatMap (\a -> [1 .. a]) (filter (>0) lambda)
  )

_convParts :: Num b => [Int] -> ([b],[b])
_convParts lambda =
  (map fromIntegral lambda, map fromIntegral (_dualPartition lambda))

hookLengths :: Fractional a => [Int] -> a -> [a]
hookLengths lambda alpha = upper ++ lower
  where
  (i,j) = _ij lambda
  (lambda', lambdaConj') = _convParts lambda
  upper = zipWith (fup lambdaConj' lambda') i j
    where
    fup x y ii jj =
      x!!(jj-1) - fromIntegral ii + alpha*(y!!(ii-1) - fromIntegral jj + 1)
  lower = zipWith (flow lambdaConj' lambda') i j
    where
    flow x y ii jj =
      x!!(jj-1) - fromIntegral ii + 1 + alpha*(y!!(ii-1) - fromIntegral jj)

_betaratio :: Fractional a => [Int] -> [Int] -> Int -> a -> a
_betaratio kappa mu k alpha = alpha * prod1 * prod2 * prod3
  where
  t = fromIntegral k - alpha * fromIntegral (mu !! (k-1))
  u = zipWith (\s kap -> t + 1 - fromIntegral s + alpha * fromIntegral kap)
               [1 .. k] (take k kappa)
  v = zipWith (\s m -> t - fromIntegral s + alpha * fromIntegral m)
               [1 .. k-1] (take (k-1) mu)
  mu' = take (mu!!(k-1)-1) (_dualPartition mu)
  w = zipWith (\s m -> fromIntegral m - t - alpha * fromIntegral s)
               [1 .. mu!!(k-1)-1] mu'
  prod1 = product $ map (\x -> x / (x+alpha-1)) u
  prod2 = product $ map (\x -> (x+alpha)/x) v
  prod3 = product $ map (\x -> (x+alpha)/x) w

_N :: [Int] -> [Int] -> Int
_N lambda mu = sum $ zipWith (*) mu prods
  where
  prods = map (\i -> product $ drop i (map (+1) lambda)) [1 .. length lambda]

jack :: forall a. (Fractional a, Ord a, IArray UArray (Maybe a))
        => [a] -> [Int] -> a -> a
jack x lambda alpha =
  case _isPartition lambda && alpha > 0 of
    False -> if _isPartition lambda
      then error "alpha must be strictly positive"
      else error "lambda is not a valid integer partition"
    True -> jac (length x) 0 lambda lambda arr0 1
      where
      nll = _N lambda lambda
      n = length x
      arr0 = listArray ((1,1), (nll, n)) (replicate (nll*n) Nothing) 
      theproduct :: Int -> a
      theproduct nu0 = if nu0 <= 1
        then 1
        else product $ map (\i -> alpha * fromIntegral i + 1) [1 .. nu0-1]
      jac :: Int -> Int -> [Int] -> [Int] -> UArray (Int,Int) (Maybe a) -> a -> a
      jac m k mu nu arr beta
        | null nu || nu!!0 == 0 || m == 0 = 1
        | length nu > m && nu!!m > 0 = 0
        | m == 1 = x!!0^(nu!!0) * theproduct (nu!!0)
        | k == 0 && isJust (arr ! (_N lambda nu, m)) =
                      fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
          s = go (jac (m-1) 0 nu nu arr 1 * beta * x!!(m-1)^(sum mu - sum nu))
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
                      go (ss + jac m ii mu nu' arr gamma) (ii+1)
                    else
                      if nu'!!0 == 0
                        then
                          go (ss + gamma * x!!(m-1)^ sum mu) (ii+1)
                        else
                          let arr' = arr // [((_N lambda nu, m), Just ss)] in
                          let jck = jac (m-1) 0 nu' nu' arr' 1 in
                          let jck' = jck * gamma *
                                       x!!(m-1)^(sum mu - sum nu') in
                          go (ss+jck') (ii+1)
                else
                  go ss (ii+1)

zonal :: (Fractional a, Ord a, IArray UArray (Maybe a)) => [a] -> [Int] -> a
zonal x lambda = c * jck
  where
  k = sum lambda
  jlambda = product (hookLengths lambda 2)
  c = 2^k * realToFrac (factorial k) / jlambda
  jck = jack x lambda 2

schur :: forall a. (Fractional a, IArray UArray (Maybe a)) => [a] -> [Int] -> a
schur x lambda =
  case _isPartition lambda of
    False -> error "lambda is not a valid integer partition"
    True -> sch n 1 lambda arr0
      where
      nll = _N lambda lambda
      n = length x
      arr0 = listArray ((1,1), (nll, length x)) (replicate (nll*n) Nothing)
      sch :: Int -> Int -> [Int] -> UArray (Int,Int) (Maybe a) -> a
      sch m k nu arr
        | null nu || nu!!0 == 0 || m == 0 = 1
        | length nu > m && nu!!m > 0 = 0
        | m == 1 = x!!0 ^ (nu!!0)
        | isJust (arr ! (_N lambda nu, m)) = fromJust $ arr ! (_N lambda nu, m)
        | otherwise = s
          where
          s = go (sch (m-1) 1 nu arr) k
          go :: Fractional a => a -> Int -> a
          go !ss ii
            | length nu < ii || nu!!(ii-1) == 0 = ss
            | otherwise =
              let u = nu!!(ii-1) in
              if length nu == ii && u > 0 || u > nu!!ii
                then
                  let nu' = (element (ii-1) .~ u-1) nu in
                  if u > 1
                    then
                      go (ss + x!!(m-1) * sch m ii nu' arr) (ii+1)
                    else
                      if nu'!!0 == 0
                        then
                           go (ss + x!!(m-1)) (ii+1)
                        else
                          let arr' = arr // [((_N lambda nu, m), Just ss)] in
                          go (ss + x!!(m-1) * sch (m-1) 1 nu' arr') (ii+1)
                else
                  go ss (ii+1)
