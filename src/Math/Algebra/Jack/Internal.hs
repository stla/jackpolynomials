{-# LANGUAGE BangPatterns #-}
module Math.Algebra.Jack.Internal
  (Partition, hookLengths, _betaratio, _isPartition, _N, skewSchurLRCoefficients)
  where
import qualified Algebra.Additive                            as AA
import qualified Algebra.Ring                                as AR
import           Data.List.Index                             ( iconcatMap )
import qualified Math.Combinat.Partitions.Integer            as MCP
import           Math.Combinat.Tableaux.LittlewoodRichardson (_lrRule)
import qualified Data.Map.Strict                             as DM

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

hookLengths :: Fractional a => Partition -> a -> [a]
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

_betaratio :: Fractional a => Partition -> Partition -> Int -> a -> a
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

(.^) :: AA.C a => Int -> a -> a
(.^) k x = if k >= 0
  then AA.sum (replicate k x)
  else AA.negate $ AA.sum (replicate (-k) x)

_fromInt :: AR.C a => Int -> a
_fromInt k = k .^ AR.one

skewSchurLRCoefficients :: AR.C a => Partition -> Partition -> DM.Map Partition a
skewSchurLRCoefficients lambda mu = 
  DM.map _fromInt $ DM.mapKeys toPartition (_lrRule lambda' mu')
  where
    toPartition :: MCP.Partition -> Partition
    toPartition (MCP.Partition part) = part 
    fromPartition :: Partition -> MCP.Partition
    fromPartition part = MCP.Partition part
    lambda' = fromPartition lambda
    mu'     = fromPartition mu