{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module JackPol
  (schurPol)
  where
import Control.Lens                             ( (.~), element )
import Data.Array                               ( Array, (!), (//), listArray )
import Data.List.Index                          (iconcatMap)
import Data.Maybe                               ( fromJust, isJust )
import Math.Combinat.Partitions.Integer.IntList (_dualPartition, _isPartition)
import MultiPol2
-- import Numeric.SpecFunctions                    (factorial)

_ij :: [Int] -> ([Int],[Int])
_ij lambda =
  (
    iconcatMap (\i a ->  replicate a (i + 1)) lambda,
    concatMap (\a -> [1 .. a]) (filter (>0) lambda)
  )

_convParts :: Num b => [Int] -> ([b],[b])
_convParts lambda =
  (map fromIntegral lambda, map fromIntegral (_dualPartition lambda))

_N :: [Int] -> [Int] -> Int
_N lambda mu = sum $ zipWith (*) mu prods
  where
  prods = map (\i -> product $ drop i (map (+1) lambda)) [1 .. length lambda]

schurPol :: Int -> [Int] -> Polynomial Int
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
