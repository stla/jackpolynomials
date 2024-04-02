{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.Jack.Internal
  (Partition
  , _productHookLengths
  , jackCoeffP
  , jackCoeffQ
  , _betaratio
  , _betaRatioOfPolynomials
  , _isPartition
  , _N
  , (.^)
  , _fromInt
  , skewSchurLRCoefficients
  , isSkewPartition)
  where
import Prelude hiding ((*), (+), (-), (/), (^), (*>), product, sum, fromIntegral)
import           Algebra.Additive                            ( (+), (-), sum )
import           Algebra.Field                               ( (/) )
import           Algebra.Ring                                ( (*), product, one )
import           Algebra.ToInteger                           ( fromIntegral )
import qualified Algebra.Additive                            as AlgAdd
import qualified Algebra.Field                               as AlgField
import qualified Algebra.Ring                                as AlgRing
import           Data.List.Index                             ( iconcatMap )
import qualified Data.Map.Strict                             as DM
import           Math.Algebra.Hspray                         ( 
                                                               RatioOfPolynomials
                                                             , Polynomial
                                                             , outerVariable
                                                             , constPoly
                                                             )
import qualified Math.Combinat.Partitions.Integer            as MCP
import           Math.Combinat.Tableaux.LittlewoodRichardson ( _lrRule )
import           Number.Ratio                                ( T( (:%) ) )

type Partition = [Int]

_isPartition :: Partition -> Bool
_isPartition []           = True
_isPartition [x]          = x > 0
_isPartition (x:xs@(y:_)) = (x >= y) && _isPartition xs

_diffSequence :: [Int] -> [Int]
_diffSequence = go where
  go (x:ys@(y:_)) = (x-y) : go ys 
  go [x]          = [x]
  go []           = []

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

_convParts :: AlgRing.C b => [Int] -> ([b], [b])
_convParts lambda =
  (map fromIntegral lambda, map fromIntegral (_dualPartition lambda))

_N :: [Int] -> [Int] -> Int
_N lambda mu = sum $ zipWith (*) mu prods
  where
  prods = map (\i -> product $ drop i (map (+1) lambda)) [1 .. length lambda]

hookLengths :: AlgField.C a => Partition -> a -> ([a], [a])
hookLengths lambda alpha = (lower, upper)
  where
    (i, j) = _ij lambda
    (lambda', lambdaConj') = _convParts lambda
    upper = zipWith (fup lambdaConj' lambda') i j
      where
        fup x y ii jj =
          x!!(jj-1) - fromIntegral ii + alpha * (y!!(ii-1) - fromIntegral (jj - 1))
    lower = zipWith (flow lambdaConj' lambda') i j
      where
        flow x y ii jj =
          x!!(jj-1) - fromIntegral ii + one + alpha * (y!!(ii-1) - fromIntegral jj)

_productHookLengths :: AlgField.C a => Partition -> a -> a
_productHookLengths lambda alpha = product lower * product upper
  where
    (lower, upper) = hookLengths lambda alpha

jackCoeffP :: AlgField.C a => Partition -> a -> a
jackCoeffP lambda alpha = one / product lower
  where
    (lower, _) = hookLengths lambda alpha

jackCoeffQ :: AlgField.C a => Partition -> a -> a
jackCoeffQ lambda alpha = one / product upper
  where
    (_, upper) = hookLengths lambda alpha

_betaratio :: AlgField.C a => Partition -> Partition -> Int -> a -> a
_betaratio kappa mu k alpha = alpha * prod1 * prod2 * prod3
  where
    mukm1 = mu !! (k-1)
    t = fromIntegral k - alpha * fromIntegral mukm1
    u = zipWith (\s kap -> t + one - fromIntegral s + alpha * fromIntegral kap)
                [1 .. k] kappa 
    v = zipWith (\s m -> t - fromIntegral s + alpha * fromIntegral m)
                [1 .. k-1] mu 
    w = zipWith (\s m -> fromIntegral m - t - alpha * fromIntegral s)
                [1 .. mukm1-1] (_dualPartition mu)
    prod1 = product $ map (\x -> x / (x + alpha - one)) u
    prod2 = product $ map (\x -> (x + alpha) / x) v
    prod3 = product $ map (\x -> (x + alpha) / x) w

_betaRatioOfPolynomials :: forall a. AlgField.C a
  => Partition -> Partition -> Int -> RatioOfPolynomials a
_betaRatioOfPolynomials kappa mu k = 
  ((x * num1 * num2 * num3) :% (den1 * den2 * den3))
  where
    mukm1 = mu !! (k-1)
    x = outerVariable :: Polynomial a
    t = constPoly (fromIntegral k) - constPoly (fromIntegral mukm1) * x
    u = zipWith 
        (
        \s kap -> 
          t - constPoly (fromIntegral $ s-1) + constPoly (fromIntegral kap) * x
        )
        [1 .. k] kappa 
    v = zipWith 
        (
        \s m -> t - constPoly (fromIntegral s) + constPoly (fromIntegral m) * x
        )
        [1 .. k-1] mu 
    w = zipWith 
        (
        \s m -> constPoly (fromIntegral m) - t - constPoly (fromIntegral s) * x
        )
        [1 .. mukm1-1] (_dualPartition mu)
    num1 = product u
    den1 = product $ map (\p -> p + x - constPoly one) u
    num2 = product $ map (\p -> p + x) v
    den2 = product v
    num3 = product $ map (\p -> p + x) w
    den3 = product w
    -- prod1 = product $ map (\x -> x / (x + alpha - one)) u
    -- prod2 = product $ map (\x -> (x + alpha) / x) v
    -- prod3 = product $ map (\x -> (x + alpha) / x) w

(.^) :: AlgAdd.C a => Int -> a -> a
(.^) k x = if k >= 0
  then AlgAdd.sum (replicate k x)
  else AlgAdd.negate $ AlgAdd.sum (replicate (-k) x)

_fromInt :: AlgRing.C a => Int -> a
_fromInt k = k .^ AlgRing.one

skewSchurLRCoefficients :: Partition -> Partition -> DM.Map Partition Int
skewSchurLRCoefficients lambda mu = 
  DM.mapKeys toPartition (_lrRule lambda' mu')
  where
    toPartition :: MCP.Partition -> Partition
    toPartition (MCP.Partition part) = part 
    fromPartition :: Partition -> MCP.Partition
    fromPartition part = MCP.Partition part
    lambda' = fromPartition lambda
    mu'     = fromPartition mu

isSkewPartition :: Partition -> Partition -> Bool
isSkewPartition lambda mu = 
  _isPartition lambda && _isPartition mu && all (>= 0) (zipWith (-) lambda mu)