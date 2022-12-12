module Math.Algebra.Jack.HypergeoPQ
  ( hypergeoPQ, _allPartitions
  ) where
import           Math.Algebra.Jack              ( zonal )

gpochhammer :: Fractional a => a -> [Int] -> a -> a
gpochhammer a kappa alpha = product $ map
  (\i -> product $ map
    (\j -> a - (fromIntegral i - 1) / alpha + fromIntegral j - 1)
    [1 .. kappa !! (i - 1)]
  )
  [1 .. length kappa]

hcoeff :: Fractional a => [a] -> [a] -> [Int] -> a -> a
hcoeff a b kappa alpha = numerator / denominator / 
  fromIntegral (factorial (sum kappa))
 where
  factorial n = product [1 .. n]
  numerator   = product $ map (\x -> gpochhammer x kappa alpha) a
  denominator = product $ map (\x -> gpochhammer x kappa alpha) b

_allPartitions :: Int -> [[Int]]
_allPartitions m = [[]] ++ (map reverse (concat ps))
 where
  ps      = [] : map parts [1 .. m]
  parts n = [n] : [ x : p | x <- [1 .. n], p <- ps !! (n - x), x <= head p ]

-- | Inefficient hypergeometric function of a matrix argument
hypergeoPQ :: (Fractional a, Ord a) => Int -> [a] -> [a] -> [a] -> a
hypergeoPQ m a b x = sum $ map (\kappa -> coeff kappa * zonal x kappa) kappas
 where
  kappas      = filter (\kap -> length kap <= length x) (_allPartitions m)
  coeff kappa = hcoeff a b kappa 2
