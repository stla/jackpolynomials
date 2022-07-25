{-# LANGUAGE BangPatterns #-}
module Internal
  where
import qualified Algebra.Ring    as AR
import           Data.List.Index ( iconcatMap )

type Partition = [Int]

_isPartition :: Partition -> Bool
_isPartition []  = True
_isPartition [x] = x > 0
_isPartition (x:xs@(y:_)) = (x >= y) && _isPartition xs

_diffSequence :: [Int] -> [Int]
_diffSequence = go where
  go (x:ys@(y:_)) = (x-y) : go ys 
  go [x] = [x]
  go []  = []

_dualPartition :: Partition -> Partition
_dualPartition [] = []
_dualPartition xs = go 0 (_diffSequence xs) [] where
  go !i (d:ds) acc = go (i+1) ds (d:acc)
  go n  []     acc = finish n acc 
  finish !j (k:ks) = replicate k j ++ finish (j-1) ks
  finish _  []     = []

_ij :: Partition -> ([Int], [Int])
_ij lambda =
  (
    iconcatMap (\i a ->  replicate a (i + 1)) lambda,
    concatMap (\a -> [1 .. a]) (filter (>0) lambda)
  )

_convParts :: Num b => [Int] -> ([b], [b])
_convParts lambda =
  (map fromIntegral lambda, map fromIntegral (_dualPartition lambda))

_N :: [Int] -> [Int] -> Int
_N lambda mu = sum $ zipWith (*) mu prods
  where
  prods = map (\i -> product $ drop i (map (+1) lambda)) [1 .. length lambda]

hookLengths :: Fractional a => [Int] -> a -> [a]
hookLengths lambda alpha = upper ++ lower
  where
    (i, j) = _ij lambda
    (lambda', lambdaConj') = _convParts lambda
    upper = zipWith (fup lambdaConj' lambda') i j
      where
        fup x y ii jj =
          x!!(jj-1) - fromIntegral ii + alpha * (y!!(ii-1) - fromIntegral jj + 1)
    lower = zipWith (flow lambdaConj' lambda') i j
      where
        flow x y ii jj =
          x!!(jj-1) - fromIntegral ii + 1 + alpha * (y!!(ii-1) - fromIntegral jj)

hookLengths' :: (Fractional a, AR.C a) => [Int] -> a -> [a]
hookLengths' = hookLengths

_betaratio :: Fractional a => [Int] -> [Int] -> Int -> a -> a
_betaratio kappa mu k alpha = alpha * prod1 * prod2 * prod3
  where
    mukm1 = mu !! (k-1)
    t = fromIntegral k - alpha * fromIntegral mukm1
    u = zipWith (\s kap -> t + 1 - fromIntegral s + alpha * fromIntegral kap)
                [1 .. k] kappa 
    v = zipWith (\s m -> t - fromIntegral s + alpha * fromIntegral m)
                [1 .. k-1] mu 
    w = zipWith (\s m -> fromIntegral m - t - alpha * fromIntegral s)
                [1 .. mukm1-1] (_dualPartition mu)
    prod1 = product $ map (\x -> x / (x + alpha - 1)) u
    prod2 = product $ map (\x -> (x + alpha) / x) v
    prod3 = product $ map (\x -> (x + alpha) / x) w

_betaratio' :: (Fractional a, AR.C a) => [Int] -> [Int] -> Int -> a -> a
_betaratio' = _betaratio
