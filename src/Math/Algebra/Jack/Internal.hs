{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Math.Algebra.Jack.Internal
  ( Partition
  , jackCoeffP
  , jackCoeffQ
  , jackCoeffC
  , jackSymbolicCoeffC
  , jackSymbolicCoeffPinv
  , jackSymbolicCoeffQinv
  , _betaratio
  , _betaRatioOfSprays
  , _isPartition
  , _N
  , _fromInt
  , skewSchurLRCoefficients
  , isSkewPartition
  , sprayToMap
  , comboToSpray
  , _kostkaNumbers
  , _inverseKostkaMatrix
  , _symbolicKostkaNumbers
  , _inverseSymbolicKostkaMatrix
  , _kostkaFoulkesPolynomial
  , _hallLittlewoodPolynomialsInSchurBasis
  , _transitionMatrixHallLittlewoodSchur
  , skewHallLittlewoodP
  , skewHallLittlewoodQ
  , flaggedSemiStandardYoungTableaux
  , tableauWeight
  , isIncreasing
  , flaggedSkewTableaux
  , skewTableauWeight
  , _skewKostkaFoulkesPolynomial
  , macdonaldPolynomialP
  , macdonaldPolynomialQ
  , skewMacdonaldPolynomialP
  , skewMacdonaldPolynomialQ
  )
  where
import           Prelude 
  hiding ((*), (+), (-), (/), (^), (*>), product, sum, fromIntegral, fromInteger, recip)
import qualified Prelude                                     as P  
import           Algebra.Additive                            ( (+), (-), sum )
import qualified Algebra.Additive                            as AlgAdd
import           Algebra.Field                               ( (/), recip )
import qualified Algebra.Field                               as AlgField
import           Algebra.Module                              ( (*>) )
import           Algebra.Ring                                ( (*), product, one
                                                             , (^), fromInteger 
                                                             )
import qualified Algebra.Ring                                as AlgRing
import           Algebra.ToInteger                           ( fromIntegral )
import qualified Data.Foldable                               as DF
import qualified Data.HashMap.Strict                         as HM
import           Data.List                                   ( 
                                                               nub
                                                             , foldl'
                                                            --  , foldl1'
                                                             , uncons
                                                             )
import           Data.List.Extra                             ( 
                                                               unsnoc
                                                             , drop1
                                                             )
import           Data.List.Index                             ( iconcatMap )
import           Data.Map.Strict                             ( Map )
import qualified Data.Map.Strict                             as DM
import           Data.Matrix                                 ( 
                                                               Matrix
                                                             , nrows
                                                             , getCol 
                                                             , getRow
                                                             , minorMatrix
                                                             , (<|>)
                                                             , (<->)
                                                             , rowVector
                                                             , colVector
                                                             , getElem
                                                             , fromLists
                                                             )
import           Data.Maybe                                  ( fromJust, isJust )
import           Data.Sequence                               ( 
                                                               Seq (..)
                                                             , (|>) 
                                                             , (<|)
                                                             , (><)
                                                             )
import qualified Data.Sequence                               as S
import           Data.Tuple.Extra                            ( fst3, both, swap )
import qualified Data.Vector                                 as V
import           Math.Algebra.Hspray                         ( 
                                                               RatioOfSprays, (%:%), (%//%), (%/%)
                                                             , unitRatioOfSprays
                                                             , zeroRatioOfSprays
                                                             , asRatioOfSprays
                                                             , Spray, (.^)
                                                             , Powers (..)
                                                             , SimpleParametricSpray
                                                             , ParametricSpray
                                                             , zeroSpray
                                                             , unitSpray
                                                             , isZeroSpray
                                                             , lone, lone'
                                                             , sumOfSprays
                                                             , productOfSprays
                                                             , FunctionLike (..)
                                                              -- , prettyParametricQSpray
                                                             )
import           Math.Combinat.Partitions.Integer            (
                                                               fromPartition
                                                             , dualPartition
                                                             , partitions
                                                             , dominates
                                                             , dominatedPartitions
                                                             , partitionWidth
                                                             , toPartitionUnsafe
                                                             , dropTailingZeros
                                                             )
-- import           Math.Combinat.Partitions.Skew               (
--                                                                mkSkewPartition
--                                                              , skewPartitionElements
--                                                              , dualSkewPartition
--                                                              , SkewPartition (..)
--                                                              )
import qualified Math.Combinat.Partitions.Integer            as MCP
-- import Math.Combinat.Partitions.Skew
-- import Math.Combinat.Partitions.Skew.Ribbon
-- import Data.Tree
import           Math.Combinat.Permutations                  ( permuteMultiset )
import           Math.Combinat.Tableaux.GelfandTsetlin       (
                                                                GT
                                                              , kostkaGelfandTsetlinPatterns
                                                             )
import           Math.Combinat.Tableaux.LittlewoodRichardson ( _lrRule )

-- f :: [Int] -> Tree MCP.Partition
-- f rho = unfoldTree (\p -> (p, map (fst . fromSkewPartition  . rbShape) $ outerRibbonsOfLength p (rho !! (partitionWidth p)))) (toPartitionUnsafe [])

-- h lambda rho=foldr (\x zs -> filter (\mu -> MCP.isSubPartitionOf mu (toPartitionUnsafe lambda)) [ribbon | z <- zs, ribbon <- map (fst . fromSkewPartition  . rbShape) $ outerRibbonsOfLength z x]) [toPartitionUnsafe []] rho
-- h' lambda rho=foldr (\x zs ->  [z++[ribbon] | z <- zs, ribbon <- map (fst . fromSkewPartition  . rbShape) $ outerRibbonsOfLength (last z) x]) [[toPartitionUnsafe []]] rho
-- -- g :: [Int] -> Tree MCP.Partition
-- -- g rho = (Tree (rootLabel (f (rho!!0))) (subForest (f (rho!!1))))

-- rhos = DM.keys $ psCombination (schurPol' (sum lambda) lambda)
-- seqs = map (\rho -> (filter ((==lambda') . last) (h' lambda rho))) rhos
-- ht ss = sum $ map (rbHeight . fromJust . toRibbon . mkSkewPartition) (zip (drop1 ss) ss)
-- e h = if even h then 1 else -1 :: Int
-- sum $ map (e . ht) seqs

-- ok q r = mu!!(q-1) - q + r > mu!!(p-1) - p
-- ok' p q r = ok q r && mu!!(p-2) - p + 1 > mu!!(q-1) - q + r
-- pairs r = map (\q -> (1, q)) (filter (\q -> ok q r) [1 .. n])
-- pairs' r = filter (\(p, q) -> ok' p q r) [(p, q) | p <- [2 .. n], q <- [2 .. n], p <= q]
-- flambda p q r = 
--   concat [
--            [mu!!(i-1) | i <- [1 .. p]]
--          , [mu!!(q-1) + p - q + r]
--          , [mu!!(i-1) + 1 | i <- [p .. q-1]]
--          , [mu!!(i-1) | i <- [q+1 .. n]]
--          ]
-- ribbons r = [ribbon (lbda \ mu) | lbda <- [flambda p q r | (p, q) <- pairs r ++ pairs' r]]
sequenceOfPartitions :: Partition -> Partition -> [[Partition]]
sequenceOfPartitions lambda rho = 
   foldr (\r zs ->  [z++[lbda ++ replicate (n-length lbda) 0] | z <- zs, lbda <- lambdas r (last z)]) [[replicate n 0]] rho
   where
    n = sum lambda
    lambda' = toPartitionUnsafe lambda
    lambdas r mu = filter (\lbda -> MCP.isSubPartitionOf (MCP.mkPartition lbda) lambda') [flambda p q r mu | (p, q) <- pairs r mu ++ pairs' r mu]
    -- filter (\lbda -> MCP.isSubPartitionOf (MCP.mkPartition lbda) lambda') 
    flambda p q r mu = 
      concat [
              [mu!!(i-1) | i <- [1 .. p-1]]
            , [mu!!(q-1) + p - q + r]
            , [mu!!(i-1) + 1 | i <- [p .. q-1]]
            , [mu!!(i-1) | i <- [q+1 .. n]]
            ]
    pairs r mu = map (\q -> (1, q)) (filter (\q -> ok q r mu) [1 .. n])
    ok q r mu = mu!!(q-1) - q + r > mu!!0 - 1
    pairs' r mu = filter (\(p, q) -> ok' p q r mu) [(p, q) | p <- [2 .. n], q <- [2 .. n], p <= q]
    ok' p q r mu = mu!!(q-1) - q + r > mu!!(p-1) - p && mu!!(p-2) - p + 1 > mu!!(q-1) - q + r

sequenceOfPartitions' :: Partition -> Partition -> [[Partition]]
sequenceOfPartitions' lambda rho = 
   foldr (\r zs ->  [z++[lbda ++ replicate (n-length lbda) 0] | z <- zs, lbda <- lambdas r (last z)]) [[replicate n 0]] rho
   where
    n = sum lambda
    lambda' = lambda ++ replicate (n - length lambda) 0
    lambdas r mu = [flambda p q r mu | (p, q) <- pairs r mu ++ pairs' r mu]
    flambda p q r mu = 
      concat [
              [mu!!(i-1) | i <- [1 .. p-1]]
            , [mu!!(q-1) + p - q + r]
            , [mu!!(i-1) + 1 | i <- [p .. q-1]]
            , [mu!!(i-1) | i <- [q+1 .. n]]
            ]
    pairs r mu = map (\q -> (1, q)) (filter (\q -> ok q r mu) [1 .. n])
    ok q r mu = mu!!(q-1) - q + r > mu!!0 - 1 && mu!!(q-1) <= lambda'!!0
    pairs' r mu = filter (\(p, q) -> ok' p q r mu) [(p, q) | p <- [2 .. n], q <- [2 .. n], p <= q]
    ok' p q r mu = mu!!(q-1) - q + r > mu!!(p-1) - p && mu!!(p-2) - p + 1 > mu!!(q-1) - q + r
                    && and (zipWith (<=) (take (p-1) mu) lambda')
                    && mu!!(q-1) <= lambda'!!(p-1)
                    && and [mu!!(i-1) + 1 <= lambda'!!i | i <- [p .. q-1]]
                    && and [mu!!(i-1) <= lambda'!!(i-1)| i <- [q+1 .. n]]

sequenceOfPartitions'' :: Seq Int -> Seq Int -> [Seq (Seq Int)]
sequenceOfPartitions'' lambda rho = 
   foldr 
     (\r zs -> 
      [z |> lbda -- (lbda >< S.replicate (n - S.length lbda) 0) 
        | z <- zs
        , lbda <- lambdas r (z `S.index` (S.length z - 1))
        , and (S.zipWith (<=) lbda lambda')
      ]) 
        [S.singleton (S.replicate n 0)] rho
   where
    n = DF.sum lambda
    lambda' = lambda >< S.replicate (n - S.length lambda) 0
    lambdas r mu = [flambda p q r mu | (p, q) <- pairs r mu ++ pairs' r mu]
    flambda p q r mu = 
      mconcat $
        map (S.fromList) 
                [
                  [mu `S.index` (i-1) | i <- [1 .. p-1]]
                , [mu `S.index` (q-1) + p - q + r]
                , [mu `S.index` (i-1) + 1 | i <- [p .. q-1]]
                , [mu `S.index` (i-1) | i <- [q+1 .. n]]
                ]
    pairs r mu = [(1, q) | q <- [1 .. n], ok q r mu]
    ok q r mu = 
      let mu_qm1 = mu `S.index` (q-1) in
        mu_qm1 - q + r > mu `S.index` 0 - 1 
          && mu_qm1 <= lambda' `S.index` 0
    pairs' r mu = 
      [(p, q) | p <- [2 .. n], q <- [p .. n], ok' p q r mu]
    ok' p q r mu = 
       let mu_qm1 = mu `S.index` (q-1) in 
        mu_qm1 - q + r > mu `S.index` (p-1) - p 
          && mu `S.index` (p-2) - p + 1 > mu_qm1 - q + r
--          && and (S.zipWith (<=) (S.take (p-1) mu) lambda')
          && mu_qm1 <= lambda' `S.index` (p-1)
          && and [mu `S.index` (i-1) + 1 <= lambda' `S.index` i | i <- [p .. q-1]]
--          && and [mu `S.index` (i-1) <= lambda' `S.index` (i-1)| i <- [q+1 .. n]]

ribbonHeight :: Seq Int -> Seq Int -> Int
ribbonHeight lambda mu = 
  DF.sum (S.zipWith (\k n -> fromEnum (k /= n)) lambda mu) + S.length lambda - S.length mu


-- pths :: Int -> Tree a -> [[a]]
-- pths n tr = go n [] tr
--   where
--     go i ps tree 
--       | i == 0 = [[rootLabel tree]]
--       | otherwise =  map (\subtree -> concat $ go (i-1) [pth ++ [rootLabel subtree] | subtree <- subForest tree, pth <- ps] subtree) (subForest tree) 

type Partition = [Int]

type PartitionsPair = (Seq Int, Seq Int)
type PairsMap = Map PartitionsPair ([(Int,Int)], [(Int,Int)]) 

-- qPochammer :: (Eq a, AlgRing.C a) => a -> a -> Int -> a
-- qPochammer q a n = 
--   AlgRing.product [AlgRing.one AlgAdd.- a AlgRing.* q AlgRing.^ (toInteger i) | i <- [0 .. n-1]]

gtPatternDiagonals' :: GT -> [Seq Int]
gtPatternDiagonals' pattern = S.empty : [diagonal j | j <- [0 .. l]]
  where
    dropTrailingZeros = S.dropWhileR (== 0)
    l = length pattern - 1
    diagonal j = 
      dropTrailingZeros
        (S.fromList
          [pattern !! r !! c | (r, c) <- zip [l-j .. l] [0 .. j]])

-- armLength :: Partition -> Map (Int, Int) Int
-- armLength lambda = DM.fromList [((i, j), lambda !! (i-1) - j) | (i, j) <- MCP.elements lambda']
--   where
--     lambda' = toPartitionUnsafe lambda

-- legLength :: Partition -> Map (Int, Int) Int
-- legLength lambda = DM.fromList [((i, j), lambda'' !! (j-1) - i) | (i, j) <- MCP.elements lambda']
--   where
--     lambda' = toPartitionUnsafe lambda
--     lambda'' = fromPartition (dualPartition lambda')

-- blambda :: (Eq a, AlgField.C a) => Partition -> (Int, Int) -> RatioOfSprays a
-- blambda lambda s = if DM.member s a then num %//% den else unitRatioOfSprays
-- blambda :: (Eq a, AlgField.C a) => Partition -> RatioOfSprays a
-- blambda lambda = num %//% den 
--   where
--     q = lone 1
--     t = lone 2
--     pairs = MCP.elements (toPartitionUnsafe lambda)
--     a = armLength lambda
--     l = legLength lambda
--     -- num = unitSpray ^-^ q^**^(a DM.! s) ^*^ t^**^(l DM.! s + 1) 
--     -- den = unitSpray ^-^ q^**^(a DM.! s + 1) ^*^ t^**^(l DM.! s) 
--     num = productOfSprays [unitSpray ^-^ q^**^(a DM.! s) ^*^ t^**^(l DM.! s + 1) | s <- pairs]
--     den = productOfSprays [unitSpray ^-^ q^**^(a DM.! s + 1) ^*^ t^**^(l DM.! s) | s <- pairs]

_dualPartition' :: Seq Int -> Seq Int
_dualPartition' Empty = S.empty
_dualPartition' xs = go 0 (_diffSequence' xs) S.empty where
  go !i (d :<| ds) acc = go (i+1) ds (d <| acc)
  go n  Empty      acc = finish n acc 
  finish !j (k :<| ks) = S.replicate k j >< finish (j-1) ks
  finish _ Empty       = S.empty
  _diffSequence' (x :<| ys@(y :<| _)) = (x-y) <| _diffSequence' ys 
  _diffSequence' (x :<| Empty)        = S.singleton x
  _diffSequence' Empty                = S.empty

-- srange :: Partition -> Partition -> [(Int, Int)]
-- srange lambda mu = [(i+1, j+1) | i <- nonEmptyRows, j <- emptyColumns] -- r \\ c
--   where
--     getSkewPartition (SkewPartition x) = x
--     skewP = mkSkewPartition (toPartitionUnsafe lambda, toPartitionUnsafe mu)
--     skewP' = dualSkewPartition skewP
--     nonEmptyRows = findIndices ((/= 0) . snd) (getSkewPartition skewP)
--     emptyColumns = findIndices ((== 0) . snd) (getSkewPartition skewP')

    -- imax = length lambda
    -- jmax = lambda !! 0
    -- elems = skewPartitionElements (mkSkewPartition (toPartitionUnsafe lambda, toPartitionUnsafe mu))
    -- is = nub $ map fst elems
    -- js = nub $ map snd elems
    -- c = [(i, j) | i <- [1 .. imax], j <- js]
    -- r = [(i, j) | i <- is, j <- [1 .. jmax]]

-- psiLambdaMu :: (Eq a, AlgField.C a) => (Partition, Partition) -> RatioOfSprays a
-- psiLambdaMu (lambda, mu) = AlgRing.product ratios
--   where
--     ss = srange lambda mu
--     ellLambda = length lambda
--     ellMu = length mu
--     mu'' = fromPartition (dualPartition (toPartitionUnsafe mu))
--     lambda'' = fromPartition (dualPartition (toPartitionUnsafe lambda))
-- --    ratio s = blambda mu s AlgField./ blambda lambda s
--     q = lone' 1
--     t = lone' 2
--     num (i, j) =
--       if i <= ellMu && j <= (mu !! (i-1))
--         then 
--           let a = mu !! (i-1) - j
--               l = mu'' !! (j-1) - i
--           in
--           (unitSpray ^-^ q a ^*^ t (l + 1))
--           %//% (unitSpray ^-^ q (a + 1) ^*^ t l)
--         else
--           unitRatioOfSprays
--     den (i, j) =
--       if i <= ellLambda && j <= (lambda !! (i-1))
--         then 
--           let a = lambda !! (i-1) - j
--               l = lambda'' !! (j-1) - i
--           in
--           (unitSpray ^-^ q a ^*^ t (l + 1))
--           %//% (unitSpray ^-^ q (a + 1) ^*^ t l)
--         else
--           unitRatioOfSprays
--     ratio s = num s AlgField./ den s
--     ratios = map ratio ss

codedRatio :: 
  PartitionsPair -> PartitionsPair -> (Int, Int) -> ([(Int,Int)], [(Int,Int)])
codedRatio (lambda, lambda') (mu, mu') (i, j)
  | i <= ellMu && j <= mu_im1 = 
      ([(a+1, l), (a', l'+1)], [(a, l+1), (a'+1, l')])
  | j <= lambda_im1 =
      ([(a', l'+1)], [(a'+1, l')])
  | otherwise = 
      ([], [])
    where
      ellMu = S.length mu
      mu_im1 = mu `S.index` (i-1)
      a = mu_im1 - j
      l = mu' `S.index` (j-1) - i
      lambda_im1 = lambda `S.index` (i-1)
      a' = lambda_im1 - j
      l' = lambda' `S.index` (j-1) - i

psiLambdaMu :: PartitionsPair -> ([(Int,Int)], [(Int,Int)])
psiLambdaMu (lambda, mu) = 
  both concat (unzip (map (swap . (codedRatio (lambda, lambda') (mu, mu'))) ss))
  where
    lambda' = _dualPartition' lambda
    mu' = _dualPartition' mu
    ellLambda = S.length lambda
    ellMu = S.length mu
    -- bools = S.zipWith (==) lambda mu >< S.replicate (S.length lambda - S.length mu) False
    -- nonEmptyRows = S.elemIndicesL False bools
    emptyRows = S.zipWith (==) lambda mu
    bools' = S.zipWith (==) lambda' mu' 
    emptyColumns = S.elemIndicesL True bools'
    ss = [
          (i+1, j+1) 
          | i <- [0 .. ellLambda - 1], i >= ellMu || not (emptyRows `S.index` i), -- not (i `S.elem` emptyRows), 
            j <- emptyColumns, j < lambda `S.index` i
        ]

-- psiLambdaMu :: (Eq a, AlgField.C a) => (Seq Int, Seq Int) -> RatioOfSprays a
-- psiLambdaMu (lambda, mu) = AlgRing.product (map ratio ss)
--   where
--     lambda' = _dualPartition' lambda
--     mu' = _dualPartition' mu
--     ellMu = S.length mu
--     bools = S.zipWith (==) lambda mu >< S.replicate (S.length lambda - ellMu) False
--     nonEmptyRows = S.elemIndicesL False bools
--     bools' = S.zipWith (==) lambda' mu' 
--     emptyColumns = S.elemIndicesL True bools'
--     ss = [(i+1, j+1) | i <- nonEmptyRows, j <- emptyColumns]
--     q = lone' 1
--     t = lone' 2
--     ratio (i, j) 
--       | i <= ellMu && j <= mu_im1 = 
--         ((unitSpray ^-^ q a ^*^ t (l + 1)) ^*^ (unitSpray ^-^ q (a' + 1) ^*^ t l'))
--           %//% ((unitSpray ^-^ q (a + 1) ^*^ t l) ^*^ (unitSpray ^-^ q a' ^*^ t (l' + 1)))
--       | j <= lambda_im1 =
--         (unitSpray ^-^ q (a' + 1) ^*^ t l')
--           %//% (unitSpray ^-^ q a' ^*^ t (l' + 1))
--       | otherwise = 
--           unitRatioOfSprays
--         where
--           mu_im1 = mu `S.index` (i-1)
--           lambda_im1 = lambda `S.index` (i-1)
--           a = mu_im1 - j
--           l = mu' `S.index` (j-1) - i
--           a' = lambda_im1 - j
--           l' = lambda' `S.index` (j-1) - i
--     -- num (i, j) =
--     --   if i <= ellMu && j <= mu `S.index` (i-1)
--     --     then 
--     --       let a = mu `S.index` (i-1) - j
--     --           l = mu'' `S.index` (j-1) - i
--     --       in
--     --       (unitSpray ^-^ q a ^*^ t (l + 1))
--     --       %//% (unitSpray ^-^ q (a + 1) ^*^ t l)
--     --     else
--     --       unitRatioOfSprays
--     -- den (i, j) =
--     --   if i <= ellLambda && j <= (lambda `S.index` (i-1))
--     --     then 
--     --       let a = lambda `S.index` (i-1) - j
--     --           l = lambda'' `S.index` (j-1) - i
--     --       in
--     --       (unitSpray ^-^ q a ^*^ t (l + 1))
--     --       %//% (unitSpray ^-^ q (a + 1) ^*^ t l)
--     --     else
--     --       unitRatioOfSprays
--     -- ratio s = num s AlgField./ den s

-- psiGT :: (Eq a, AlgField.C a) => GT -> RatioOfSprays a
-- psiGT pattern = 
--   AlgRing.product rOS
-- --  productOfSprays numFactors %//% productOfSprays denFactors
--   where
--     lambdas = gtPatternDiagonals' pattern
--     pairs = zip (drop1 lambdas) lambdas
--     rOS = map psiLambdaMu pairs
-- --    rOS = map psiLambdaMu pairs
--     -- ell = length lambdas - 1
--     -- q = lone 1
--     -- t = lone 2
--     -- poch = qPochammer q
--     -- pairs = zip lambdas (drop1 lambdas)
--     -- (numFactors, denFactors) = unzip [
--     --     let (lambda1, lambda2) = pairs !! (k-1)
--     --         l1i = lambda1 !! i
--     --         l1j = lambda1 !! j
--     --         l2i = lambda2 !! i
--     --         l2jp1 = lambda2 !! (j + 1)
--     --     in
--     --       (
--     --         poch (q^**^(l1i - l1j) ^*^ t^**^(j-i+1)) (l2i - l1i) 
--     --         ^*^ poch (q^**^(l1i - l2jp1 + 1) ^*^ t^**^(j-i)) (l2i - l1i)
--     --       , poch (q^**^(l1i - l1j + 1) ^*^ t^**^(j-i)) (l2i - l1i) 
--     --         ^*^ poch (q^**^(l1i - l2jp1) ^*^ t^**^(j-i+1)) (l2i - l1i)
--     --       )
--     --     | k <- [1 .. ell], j <- [0 .. k-1], i <- [0 .. j]
--     --   ]
--     -- -- numFactors = [
--     -- --     let lambda1 = lambdas !! (k-1)
--     -- --         lambda2 = lambdas !! k
--     -- --     in
--     -- --       poch (q^**^(lambda1 !! i - lambda1 !! j) ^*^ t^**^(j-i+1)) (lambda2 !! i - lambda1 !! i) 
--     -- --       ^*^ poch (q^**^(lambda1 !! i - lambda2 !! (j+1) + 1) ^*^ t^**^(j-i)) (lambda2 !! i - lambda1 !! i)
--     -- --     | k <- [1 .. ell], j <- [0 .. k-1], i <- [0 .. j]
--     -- --   ]
--     -- -- denFactors = [
--     -- --     let lambda1 = lambdas !! (k-1)
--     -- --         lambda2 = lambdas !! k
--     -- --     in
--     -- --     | k <- [1 .. ell], j <- [0 .. k-1], i <- [0 .. j]
--     -- --   ]
-- -- \psi_T = \prod_{1 \le i \le j \le k - 1 \le \ell - 1} 
-- -- {(q^{\lambda^{k - 1}_i - \lambda^{k - 1}_j} t^{j - i + 1})_{\lambda^k_i - \lambda^{k - 1}_i} 
-- --  (q^{\lambda^{k - 1}_i - \lambda^k_{j + 1} + 1} t^{j - i})_{\lambda^k_i - \lambda^{k - 1}_i} 
-- -- \over (q^{\lambda^{k - 1}_i - \lambda^{k - 1}_j + 1} t^{j - i})_{\lambda^k_i - \lambda^{k - 1}_i} 
-- -- (q^{\lambda^{k - 1}_i - \lambda^k_{j + 1}} t^{j - i + 1})_{\lambda^k_i - \lambda^{k - 1}_i}}.

-- phiLambdaMu' :: (Eq a, AlgField.C a) => (Seq Int, Seq Int) -> RatioOfSprays a
-- phiLambdaMu' (lambda, mu) = bb AlgRing.* psiLambdaMu' (lambda, mu)
--   where
--     bb = (blambda (DF.toList lambda)) AlgField./ (blambda (DF.toList mu))

-- phiLambdaMu :: (Eq a, AlgField.C a) => (Seq Int, Seq Int) -> RatioOfSprays a
-- phiLambdaMu (lambda, mu) = AlgRing.product (map ratio ss)
--   where
--     lambda' = _dualPartition' lambda
--     mu' = _dualPartition' mu
--     bools' = 
--       S.zipWith (==) lambda' mu' 
--         >< S.replicate (S.length lambda' - S.length mu') False 
--     nonEmptyColumns = S.elemIndicesL False bools'
--     ss = [(i, j+1) | j <- nonEmptyColumns, i <- [1 .. lambda' `S.index` j]] 
--     q = lone' 1
--     t = lone' 2
--     ellMu = S.length mu
--     ratio (i, j) 
--       | i <= ellMu && j <= mu_im1 = 
--         ((unitSpray ^-^ q (a + 1) ^*^ t l) 
--             ^*^ (unitSpray ^-^ q a' ^*^ t (l' + 1)))
--               %//% ((unitSpray ^-^ q a ^*^ t (l + 1)) 
--                       ^*^ (unitSpray ^-^ q (a' + 1) ^*^ t l'))
--       | j <= lambda_im1 =
--           (unitSpray ^-^ q a' ^*^ t (l' + 1))
--             %//% (unitSpray ^-^ q (a' + 1) ^*^ t l') 
--       | otherwise = 
--           unitRatioOfSprays
--         where
--           mu_im1 = mu `S.index` (i-1)
--           a = mu_im1 - j
--           l = mu' `S.index` (j-1) - i
--           lambda_im1 = lambda `S.index` (i-1)
--           a' = lambda_im1 - j
--           l' = lambda' `S.index` (j-1) - i

phiLambdaMu :: PartitionsPair -> ([(Int,Int)], [(Int,Int)])
phiLambdaMu (lambda, mu) = 
  both concat (unzip (map (codedRatio (lambda, lambda') (mu, mu')) ss))
  where
    lambda' = _dualPartition' lambda
    mu' = _dualPartition' mu
    bools' = 
      S.zipWith (==) lambda' mu' 
        >< S.replicate (S.length lambda' - S.length mu') False 
    nonEmptyColumns = S.elemIndicesL False bools'
    ss = [(i, j+1) | j <- nonEmptyColumns, i <- [1 .. lambda' `S.index` j]] 
    -- ellMu = S.length mu
    -- ratio (i, j) 
    --   | i <= ellMu && j <= mu_im1 = 
    --       ([(a+1, l), (a', l'+1)], [(a,l+1), (a'+1, l)])
    --   | j <= lambda_im1 =
    --       ([(a', l'+1)], [(a'+1, l')])
    --   | otherwise = 
    --       ([], [])
    --     where
    --       mu_im1 = mu `S.index` (i-1)
    --       a = mu_im1 - j
    --       l = mu' `S.index` (j-1) - i
    --       lambda_im1 = lambda `S.index` (i-1)
    --       a' = lambda_im1 - j
    --       l' = lambda' `S.index` (j-1) - i

-- phiGT''' :: (Eq a, AlgField.C a) => GT -> RatioOfSprays a
-- phiGT''' pattern = num %//% den
--   where
--     lambdas = gtPatternDiagonals' pattern
--     pairs = zip (drop1 lambdas) lambdas
--     (num_als, den_als) = both concat (unzip (map phiLambdaMu' pairs))
--     num_als' = num_als \\ den_als
--     den_als' = den_als \\ num_als
--     num_assocs = DM.assocs (foldl' (\m k -> DM.insertWith (+) k 1 m) DM.empty num_als')
--     den_assocs = DM.assocs (foldl' (\m k -> DM.insertWith (+) k 1 m) DM.empty den_als')
--     q = lone' 1
--     t = lone' 2
--     poly ((a, l), c) = (unitSpray ^-^ q a ^*^ t l) ^**^ c
--     num = productOfSprays (map poly num_assocs)
--     den = productOfSprays (map poly den_assocs)

makeRatioOfSprays :: 
  (Eq a, AlgField.C a) => 
  PairsMap -> [PartitionsPair] -> RatioOfSprays a
makeRatioOfSprays pairsMap pairs = num %//% den
  where
    als = both concat (unzip (map ((DM.!) pairsMap) pairs))
    (num_map, den_map) =
      both (foldl' (\m k -> DM.insertWith (+) k 1 m) DM.empty) als
    f k1 k2 = if k1 > k2 then Just (k1 - k2) else Nothing
    assocs = both DM.assocs
      (
        DM.differenceWith f num_map den_map
      , DM.differenceWith f den_map num_map
      )
    -- (num_als, den_als) = both concat (unzip (map ((DM.!) pairsMap) pairs))
    -- als = (num_als \\ den_als, den_als \\ num_als)
    -- assocs = 
    --   both (DM.assocs . (foldl' (\m k -> DM.insertWith (+) k 1 m) DM.empty)) als
    q = lone' 1
    t = lone' 2
    poly ((a, l), c) = (unitSpray ^-^ q a ^*^ t l) ^**^ c
    (num, den) = both (productOfSprays . (map poly)) assocs

_macdonaldPolynomial :: 
  (Eq a, AlgField.C a) 
  => (PartitionsPair -> ([(Int,Int)], [(Int,Int)]))
  -> Int 
  -> Partition 
  -> ParametricSpray a
_macdonaldPolynomial f n lambda = HM.unions hashMaps
  where
    lambda' = toPartitionUnsafe lambda
    mus = filter (\mu -> partitionWidth mu <= n) (dominatedPartitions lambda')
    pairing lambdas = zip (drop1 lambdas) lambdas
    listsOfPairs = 
      map (
        map (pairing . gtPatternDiagonals') 
          . (kostkaGelfandTsetlinPatterns lambda')
      ) mus
    allPairs = nub $ concat (concat listsOfPairs)
    pairsMap = DM.fromList (zip allPairs (map f allPairs))
    coeffs = HM.fromList 
      (zipWith 
        (\mu listOfPairs -> 
          (
            S.fromList (fromPartition mu)
          , AlgAdd.sum (map (makeRatioOfSprays pairsMap) listOfPairs)
          )
        ) mus listsOfPairs
      )
    dropTrailingZeros = S.dropWhileR (== 0)
    hashMaps = 
      map 
        (\mu -> 
          let mu' = fromPartition mu
              mu'' = S.fromList mu'
              mu''' = mu' ++ (replicate (n - S.length mu'') 0)
              coeff = coeffs HM.! mu''
              compos = permuteMultiset mu'''
          in
            HM.fromList 
              [let compo' = dropTrailingZeros (S.fromList compo) in
                (Powers compo' (S.length compo'), coeff) | compo <- compos]
        ) mus


-- phiGT :: (Eq a, AlgField.C a) => GT -> RatioOfSprays a
-- phiGT pattern = 
--   AlgRing.product rOS
--   where
--     lambdas = gtPatternDiagonals' pattern
--     pairs = zip (drop1 lambdas) lambdas
--     rOS = map phiLambdaMu pairs
  -- productOfSprays numFactors %//% productOfSprays denFactors
  -- where
  --   lambdas = gtPatternDiagonals' pattern
  --   ell = length lambdas - 1
  --   q = lone 1
  --   t = lone 2
  --   poch = qPochammer q
  --   numFactors = [
  --       let lambda1 = lambdas !! (k-1)
  --           lambda2 = lambdas !! k
  --       in
  --         poch (q^**^(lambda2 !! i - lambda2 !! j) ^*^ t^**^(j-i+1)) (lambda2 !! j - lambda1 !! j) 
  --         ^*^ poch (q^**^(lambda1 !! i - lambda2 !! (j+1) + 1) ^*^ t^**^(j-i)) (lambda2 !! j - lambda1 !! j)
  --       | k <- [1 .. ell], j <- [0 .. k], i <- [0 .. j]
  --     ]
  --   denFactors = [
  --       let lambda1 = lambdas !! (k-1)
  --           lambda2 = lambdas !! k
  --       in
  --         poch (q^**^(lambda2 !! i - lambda2 !! j + 1) ^*^ t^**^(j-i)) (lambda2 !! j - lambda1 !! j) 
  --         ^*^ poch (q^**^(lambda2 !! i - lambda2 !! j + 1) ^*^ t^**^(j-i)) (lambda2 !! j - lambda1 !! j)
  --       | k <- [1 .. ell], j <- [0 .. k], i <- [0 .. j]
  --     ]
-- \phi_T = \prod_{1 \le i \le j \le k \le \ell} 
--{(q^{\lambda^k_i - \lambda^k_j} t^{j - i + 1})_{\lambda^k_j - \lambda^{k - 1}_j} 
-- (q^{\lambda^k_i - \lambda^k_{j + 1} + 1} t^{j - i})_{\lambda^k_{j + 1} - \lambda^{k - 1}_{j + 1}} 
-- \over (q^{\lambda^k_i - \lambda^k_j + 1} t^{j - i})_{\lambda^k_j - \lambda^{k - 1}_j} 
-- (q^{\lambda^k_i - \lambda^k_{j + 1}} t^{j - i + 1})_{\lambda^k_{j + 1} -\lambda^{k - 1}_{j + 1}}}.

-- macdonaldPolynomialP :: 
--   (Eq a, AlgField.C a) => Int -> Partition -> ParametricSpray a
-- macdonaldPolynomialP n lambda = HM.unions hashMaps -- sumOfSprays sprays
--   where
--     lambda' = toPartitionUnsafe lambda
--     mus = filter (\mu -> partitionWidth mu <= n) (dominatedPartitions lambda')
--     dropTrailingZeros = S.dropWhileR (== 0)
--     coeffs = HM.fromList 
--       [
--         (
--           S.fromList (fromPartition mu)
--         , AlgAdd.sum (map psiGT (kostkaGelfandTsetlinPatterns lambda' mu))
--         ) | mu <- mus
--       ]
--     hashMaps = 
--       map 
--         (\mu -> 
--           let mu' = fromPartition mu
--               mu'' = S.fromList mu'
--               mu''' = mu' ++ (replicate (n - S.length mu'') 0)
--               coeff = coeffs HM.! mu''
--               compos = permuteMultiset mu'''
--           in
--             HM.fromList 
--               [let compo' = dropTrailingZeros (S.fromList compo) in
--                 (Powers compo' (S.length compo'), coeff) | compo <- compos]
--         ) mus
--     -- compos = compositions n (sum lambda)
--     -- compos1 = map (filter (/= 0)) compos
--     -- lambda' = toPartitionUnsafe lambda
--     -- coeffs = DM.fromList 
--     --   [
--     --     (
--     --       compo1
--     --     , AlgAdd.sum (map psiGT (kostkaGelfandTsetlinPatterns' lambda' compo1))
--     --     ) | compo1 <- nub compos1
--     --   ]
--     -- zippedCompos = zip compos compos1
--     -- term (compo, compo1) = 
--     --   coeffs DM.! compo1 *^ monomial (filter ((/= 0) . snd) $ zip [1 ..] compo)
--     -- sprays = map term zippedCompos

-- test :: Bool
-- test = 
--   prettyParametricQSpray (HM.map (\rOQ -> substitute [Just 0, Nothing] rOQ) macPoly)
--     == "{ [ 1 ] }*X^2.Y + { [ 1 ] }*X^2.Z + { [ 1 ] }*X.Y^2 + { [ -a2^2 - a2 + 2 ] }*X.Y.Z + { [ 1 ] }*X.Z^2 + { [ 1 ] }*Y^2.Z + { [ 1 ] }*Y.Z^2"
--   where
--     macPoly = macdonaldPolynomialP 3 [2, 1]

macdonaldPolynomialP :: 
  (Eq a, AlgField.C a) => Int -> Partition -> ParametricSpray a
macdonaldPolynomialP = _macdonaldPolynomial psiLambdaMu 

macdonaldPolynomialQ :: 
  (Eq a, AlgField.C a) => Int -> Partition -> ParametricSpray a
macdonaldPolynomialQ = _macdonaldPolynomial phiLambdaMu 
-- macdonaldPolynomialQ n lambda = HM.unions hashMaps
--   where
--     lambda' = toPartitionUnsafe lambda
--     mus = filter (\mu -> partitionWidth mu <= n) (dominatedPartitions lambda')
--     pairing lambdas = zip (drop1 lambdas) lambdas
--     listsOfPairs = 
--       map (
--         map (pairing . gtPatternDiagonals') 
--           . (kostkaGelfandTsetlinPatterns lambda')
--       ) mus
--     allPairs = nub $ concat (concat listsOfPairs)
--     pairsMap = DM.fromList (zip allPairs (map phiLambdaMu allPairs))
--     coeffs = HM.fromList 
--       (zipWith 
--         (\mu listOfPairs -> 
--           (
--             S.fromList (fromPartition mu)
--           , AlgAdd.sum (map (phiGT pairsMap) listOfPairs)
--           )
--         ) mus listsOfPairs
--       )
--     dropTrailingZeros = S.dropWhileR (== 0)
--     hashMaps = 
--       map 
--         (\mu -> 
--           let mu' = fromPartition mu
--               mu'' = S.fromList mu'
--               mu''' = mu' ++ (replicate (n - S.length mu'') 0)
--               coeff = coeffs HM.! mu''
--               compos = permuteMultiset mu'''
--           in
--             HM.fromList 
--               [let compo' = dropTrailingZeros (S.fromList compo) in
--                 (Powers compo' (S.length compo'), coeff) | compo <- compos]
--         ) mus

-- phiGT' :: (Eq a, AlgField.C a) => Map (Seq Int, Seq Int) (RatioOfSprays a) -> GT -> RatioOfSprays a
-- phiGT' pairsMap pattern = 
--   AlgRing.product rOS
--   where
--     lambdas = gtPatternDiagonals' pattern
--     pairs = zip (drop1 lambdas) lambdas
--     rOS = map ((DM.!) pairsMap) pairs

-- phiGT'' :: 
--   (Eq a, AlgField.C a) => 
--   Map (Seq Int, Seq Int) (RatioOfSprays a) -> [(Seq Int, Seq Int)] -> RatioOfSprays a
-- phiGT'' pairsMap pairs = 
--   AlgRing.product rOS
--   where
--     rOS = map ((DM.!) pairsMap) pairs

-- tt :: Int -> Partition -> [[[(Seq Int, Seq Int)]]]
-- tt n lambda = listsOfPairs
--   where
--     lambda' = toPartitionUnsafe lambda
--     mus = filter (\mu -> partitionWidth mu <= n) (dominatedPartitions lambda')
--     pairs lambdas = zip (drop1 lambdas) lambdas
--     listsOfPairs = 
--       map (
--         map (pairs . gtPatternDiagonals') 
--           . (kostkaGelfandTsetlinPatterns lambda')
--       ) mus

-- macdonaldPolynomialQ :: 
--   (Eq a, AlgField.C a) => Int -> Partition -> ParametricSpray a
-- macdonaldPolynomialQ n lambda = HM.unions hashMaps -- sumOfSprays sprays
--   where
--     lambda' = toPartitionUnsafe lambda
--     mus = filter (\mu -> partitionWidth mu <= n) (dominatedPartitions lambda')
--     pairing lambdas = zip (drop1 lambdas) lambdas
--     listsOfPairs = 
--       map (
--         map (pairing . gtPatternDiagonals') 
--           . (kostkaGelfandTsetlinPatterns lambda')
--       ) mus
--     -- allPairs = map (\mu -> map (pairs . gtPatternDiagonals') (kostkaGelfandTsetlinPatterns lambda' mu)) mus
-- --      map (\mu -> map (pairs . gtPatternDiagonals') (kostkaGelfandTsetlinPatterns lambda' mu)) mus
--     -- allPairs =
--     --   map 
--     --     (concatMap (pairs . gtPatternDiagonals') . (kostkaGelfandTsetlinPatterns lambda')) mus
--     -- allPairs =
--     --   concatMap 
--     --     (concatMap (pairs . gtPatternDiagonals') . (kostkaGelfandTsetlinPatterns lambda')) mus
--     allPairs = nub $ concat (concat listsOfPairs)
--     pairsMap = DM.fromList (zip allPairs (map phiLambdaMu allPairs))
--     coeffs = HM.fromList 
--       (zipWith 
--         (\mu listOfPairs -> 
--           (
--             S.fromList (fromPartition mu)
--           , AlgAdd.sum (map (phiGT'' pairsMap) listOfPairs)
--           )
--         ) mus listsOfPairs
--       )
--     -- allPairs = nub $ 
--     --   concatMap 
--     --     (concatMap (pairs . gtPatternDiagonals') . (kostkaGelfandTsetlinPatterns lambda')) mus
--     -- pairsMap = DM.fromList (zip allPairs (map phiLambdaMu allPairs))
--     -- coeffs = HM.fromList 
--     --   [
--     --     (
--     --       S.fromList (fromPartition mu)
--     --     , AlgAdd.sum (map (phiGT' pairsMap) (kostkaGelfandTsetlinPatterns lambda' mu))
--     --     ) | mu <- mus
--     --   ]
--     dropTrailingZeros = S.dropWhileR (== 0)
--     hashMaps = 
--       map 
--         (\mu -> 
--           let mu' = fromPartition mu
--               mu'' = S.fromList mu'
--               mu''' = mu' ++ (replicate (n - S.length mu'') 0)
--               coeff = coeffs HM.! mu''
--               compos = permuteMultiset mu'''
--           in
--             HM.fromList 
--               [let compo' = dropTrailingZeros (S.fromList compo) in
--                 (Powers compo' (S.length compo'), coeff) | compo <- compos]
--         ) mus

-- test' :: Bool 
-- test' = 
--   prettyParametricQSpray (HM.map (\rOQ -> substitute [Just 0, Nothing] rOQ) macPoly)
--     == "{ [ a2^2 - 2*a2 + 1 ] }*X^2.Y + { [ a2^2 - 2*a2 + 1 ] }*X^2.Z + { [ a2^2 - 2*a2 + 1 ] }*X.Y^2 + { [ -a2^4 + a2^3 + 3*a2^2 - 5*a2 + 2 ] }*X.Y.Z + { [ a2^2 - 2*a2 + 1 ] }*X.Z^2 + { [ a2^2 - 2*a2 + 1 ] }*Y^2.Z + { [ a2^2 - 2*a2 + 1 ] }*Y.Z^2"
--   where
--     macPoly = macdonaldPolynomialQ 3 [2, 1]

lastSubPartition :: Int -> Partition -> Partition
-- assumes w <= sum(k:ks)
lastSubPartition 0 _  = []
lastSubPartition _ [] = []
lastSubPartition w (k:ks) =  
  if w <= k then [w] else k : lastSubPartition (w - k) ks

_skewMacdonaldPolynomial :: 
  (Eq a, AlgField.C a) 
  => (PartitionsPair -> ([(Int,Int)], [(Int,Int)]))
  -> Int 
  -> Partition 
  -> Partition
  -> ParametricSpray a
_skewMacdonaldPolynomial f n lambda mu = HM.unions hashMaps
  where
    nus = 
      filter ((<= n) . partitionWidth) $ 
        dominatedPartitions 
          (toPartitionUnsafe (lastSubPartition (sum lambda - sum mu) lambda))
    pairing lambdas = zip (drop1 lambdas) lambdas
    mapOfPatterns = HM.filter (not . null) 
      (HM.fromList (map (\nu -> 
        let nu' = fromPartition nu in
          (
            S.fromList nu'
          , skewGelfandTsetlinPatterns lambda mu nu'
          )        
        ) nus))
    mapOfPairs = HM.map (map pairing) mapOfPatterns
    listsOfPairs = HM.elems mapOfPairs
    allPairs = nub $ concat (concat listsOfPairs)
    pairsMap = DM.fromList (zip allPairs (map f allPairs))
    coeffs = 
      HM.map (AlgAdd.sum . (map (makeRatioOfSprays pairsMap))) mapOfPairs
    dropTrailingZeros = S.dropWhileR (== 0)
    hashMaps = 
      map 
        (\nu'' -> 
          let nu''' = DF.toList nu'' ++ (replicate (n - S.length nu'') 0)
              coeff = coeffs HM.! nu''
              compos = permuteMultiset nu'''
          in
            HM.fromList 
              [let compo' = dropTrailingZeros (S.fromList compo) in
                (Powers compo' (S.length compo'), coeff) | compo <- compos]
        ) (HM.keys coeffs)

-- xxx :: Int -> Partition -> Partition -> (RatioOfSprays Rational, RatioOfSprays Rational)
-- xxx n lambda mu = (roq1, roq2)
--   where
--     nus = partitions' (lambda !! 0, n) (sum lambda - sum mu)
--     pairing lambdas = zip (drop1 lambdas) lambdas
--     listsOfPairs = 
--       map (
--         map pairing 
--           . (skewGelfandTsetlinPatterns lambda mu)
--           . fromPartition
--       ) nus
--     allPairs = nub $ concat (concat listsOfPairs)
--     pairsMap = DM.fromList (zip allPairs (map psiLambdaMu allPairs))
--     roq1 = makeRatioOfSprays pairsMap (listsOfPairs !! 0 !! 0)
--     roq2 = makeRatioOfSprays pairsMap (listsOfPairs !! 0 !! 1)

skewMacdonaldPolynomialP :: 
  (Eq a, AlgField.C a) => Int -> Partition -> Partition -> ParametricSpray a
skewMacdonaldPolynomialP = _skewMacdonaldPolynomial psiLambdaMu 

skewMacdonaldPolynomialQ :: 
  (Eq a, AlgField.C a) => Int -> Partition -> Partition -> ParametricSpray a
skewMacdonaldPolynomialQ = _skewMacdonaldPolynomial phiLambdaMu 

-- test :: String
-- test = 
--   prettyParametricQSpray (HM.map (\rOQ -> substitute [Just 0, Nothing] rOQ) macPoly)
-- --    == "{ [ 1 ] }*X^2.Y + { [ 1 ] }*X^2.Z + { [ 1 ] }*X.Y^2 + { [ -a2^2 - a2 + 2 ] }*X.Y.Z + { [ 1 ] }*X.Z^2 + { [ 1 ] }*Y^2.Z + { [ 1 ] }*Y.Z^2"
--   where
--     macPoly = skewMacdonaldPolynomialP 3 [3, 2] [1, 1]

sandwichedPartitions :: Int -> Seq Int -> Seq Int -> [Seq Int]
sandwichedPartitions weight mu lambda = 
  recursiveFun weight (lambda `S.index` 0) mu' lambda
  where
    mu' = mu >< (S.replicate (S.length lambda - S.length mu) 0)
    dropTrailingZeros = S.dropWhileR (== 0)
    recursiveFun :: Int -> Int -> Seq Int -> Seq Int -> [Seq Int]
    recursiveFun d h0 a_as b_bs
      | d < 0 || d < DF.sum a_as || d > DF.sum b_bs = []
      | d == 0 = [S.empty]
      | otherwise = 
          concatMap 
            (\h -> 
              [h :<| dropTrailingZeros hs | hs <- recursiveFun (d-h) h as bs]
            )
            [max 0 a .. min h0 b]
          where
            a = a_as `S.index` 0
            b = b_bs `S.index` 0
            as = S.drop 1 a_as
            bs = S.drop 1 b_bs

skewGelfandTsetlinPatterns :: Partition -> Partition -> [Int] -> [[Seq Int]]
skewGelfandTsetlinPatterns lambda mu weight 
  -- | not (isSkewPartition lambda mu) =
  --     error "skewGelfandTsetlinPatterns: invalid skew partition."
  | any (< 0) weight =
      []
  | wWeight /= wLambda - wMu = 
      []
  | wWeight == 0 =
      [replicate (length weight + 1) lambda']
  | otherwise =
      if any (== 0) weight 
        then map (\pattern -> [pattern `S.index` i | i <- indices]) patterns
        else map DF.toList patterns
  where
    wWeight = sum weight
    lambda' = S.fromList lambda
    wLambda = DF.sum lambda'
    mu' = S.fromList mu
    wMu = DF.sum mu'
    recursiveFun :: Seq Int -> Seq Int -> [Seq (Seq Int)]
    recursiveFun kappa w =
      if d == wMu 
        then
          if ellKappa >= ellMu && 
              and (S.zipWith (>=) kappa mu') && 
                ellKappa <= ellMu + 1 && 
                  and (S.zipWith (>=) (mu') (S.drop 1 kappa))
            then [S.fromList [mu', kappa]]
            else [] 
        else 
          concatMap
            (\nu -> [list |> kappa | list <- recursiveFun nu hw])
              (sandwichedPartitions d (S.drop 1 kappa |> 0) kappa)
        where
          ellKappa = S.length kappa
          ellMu = S.length mu'
          d = DF.sum kappa - w `S.index` 0
          hw = S.drop 1 w
    weight' = S.filter (/= 0) (S.fromList weight)
    patterns = recursiveFun lambda' (S.reverse weight')
    indices = map (subtract 1) (scanl1 (+) (1 : map (min 1) (reverse weight)))

skewGelfandTsetlinPatternToTableau :: [Seq Int] -> [(Int, Seq Int)]
skewGelfandTsetlinPatternToTableau pattern = 
  if ellLambda == 0
    then []
    else DF.toList skewTableau
  where
    lambda = pattern !! (length pattern - 1)
    ellLambda = S.length lambda
    mu = pattern !! 0
    mu' = mu >< (S.replicate (ellLambda - S.length mu) 0)
    skewPartitionRows kappa nu = 
      concatMap (uncurry replicate) (S.zip differences indices)
      where
        indices = S.fromList [0 .. ellLambda]
        differences = S.zipWith (-) kappa nu >< S.drop (S.length nu) kappa
    startingTableau = S.replicate ellLambda S.Empty
    growTableau :: Seq (Seq Int) -> (Int, Seq Int, Seq Int) -> Seq (Seq Int)
    growTableau tableau (j, kappa, nu) =
      DF.foldr (S.adjust' (flip (|>) j)) tableau (skewPartitionRows kappa nu)
    skewPartitions = zip3 [1 ..] (drop1 pattern) pattern
    skewTableau = 
      S.zip mu' (DF.foldl' growTableau startingTableau skewPartitions)

skewTableauxWithGivenShapeAndWeight :: 
  Partition -> Partition -> [Int] -> [[(Int, Seq Int)]]
skewTableauxWithGivenShapeAndWeight lambda mu weight = 
  map skewGelfandTsetlinPatternToTableau 
      (skewGelfandTsetlinPatterns lambda mu weight) 

_skewKostkaFoulkesPolynomial :: 
  (Eq a, AlgRing.C a) => Partition -> Partition -> Partition -> Spray a
_skewKostkaFoulkesPolynomial lambda mu nu = 
  if sum lambda == sum mu + sum nu
    then sumOfSprays sprays
    else zeroSpray
  where
    tableaux = skewTableauxWithGivenShapeAndWeight lambda mu nu
    word skewT = mconcat (map S.reverse (snd (unzip skewT))) 
    mm = lone' 1 
    sprays = map (mm . charge . word) tableaux

gtPatternDiagonals :: GT -> (Int, [Partition])
gtPatternDiagonals pattern = (corner, [diagonal j | j <- [1 .. l]])
  where
    l = length pattern - 1
    corner = pattern !! l !! 0
    diagonal j = 
      dropTailingZeros
        [pattern !! r !! c | (r, c) <- zip [l-j .. l] [0 .. j]]

gtPatternToTableau :: GT -> [Seq Int]
gtPatternToTableau pattern = 
  if l >= 0 
    then DF.toList $ go 0 startingTableau
    else [S.replicate corner 1]
  where
    (corner, diagonals) = gtPatternDiagonals pattern
    diagonals' = [corner] : diagonals
    l = length diagonals - 1
    lambda = diagonals !! l
    m = length lambda
    startingTableau = S.replicate m S.Empty
    skewPartitions = zip diagonals diagonals'
    skewPartitionRows (kappa, nu) = 
      concatMap (\(i, d) -> replicate d i) (zip [0 ..] differences)
      where
        differences = zipWith (-) kappa nu ++ drop (length nu) kappa
    go i tableau
      | i == 0 = go 1 (S.adjust' (flip (><) (S.replicate corner 1)) 0 tableau)
      | i == l+2 = tableau
      | otherwise = 
          go (i+1) (growTableau (i+1) tableau (skewPartitions !! (i-1)))
    growTableau :: 
      Int -> Seq (Seq Int) -> (Partition, Partition) -> Seq (Seq Int)
    growTableau j tableau skewPart =
      DF.foldr (S.adjust' (flip (|>) j)) tableau (skewPartitionRows skewPart)

semiStandardTableauxWithGivenShapeAndWeight :: 
  Partition -> Partition -> [[Seq Int]]
semiStandardTableauxWithGivenShapeAndWeight lambda mu =
  if lambda' `dominates` mu'
    then map gtPatternToTableau (kostkaGelfandTsetlinPatterns lambda' mu')
    else []
  where
    lambda' = toPartitionUnsafe lambda
    mu' = toPartitionUnsafe mu

-- length lambda = length as = length bs; as <= bs; last bs >= length lambda
flaggedSemiStandardYoungTableaux :: Partition -> [Int] -> [Int] -> [[[Int]]] 
flaggedSemiStandardYoungTableaux lambda as bs = 
  worker (repeat 0) lambda 0
    where
      worker _ [] _ = [[]] 
      worker prevRow (s:ss) i
        = [ (r:rs) 
            | r <- row (bs !! i) s (as !! i) prevRow
            , rs <- worker (map (+1) r) ss (i + 1) ]
      -- weekly increasing lists of length @len@, pointwise at least @xs@, 
      -- maximum value @n@, minimum value @prev@.
      row :: Int -> Int -> Int -> [Int] -> [[Int]]
      row n len prev xxs = 
        if len == 0 
          then [[]] 
          else [ (j:js) | j <- [max x prev .. n], js <- row n (len-1) j xs ]
          where
            (x, xs) = fromJust (uncons xxs)

tableauWeight :: [[Int]] -> [Int]
tableauWeight tableau = [count i | i <- [1 .. m]]
  where
    x = concat tableau
    m = maximum x
    count i = sum [fromEnum (k == i) | k <- x]

flaggedSkewTableaux :: 
  Partition -> Partition -> [Int] -> [Int] -> [[(Int,[Int])]]
flaggedSkewTableaux lambda mu as bs = worker uus vvs dds (repeat 1) 0
  where
    uus = mu ++ (replicate (length lambda - length mu) 0)
    vvs = zipWith (-) lambda uus
    dds = _diffSequence uus
    _diffSequence :: [Int] -> [Int]
    _diffSequence = go where
      go (x:ys@(y:_)) = (x-y) : go ys 
      go [x] = [x]
      go []  = []
    -- | @worker inner outerMinusInner innerdiffs lowerbound
    worker :: [Int] -> [Int] -> [Int] -> [Int] -> Int -> [[(Int,[Int])]]
    worker (u:us) (v:vs) (d:ds) lb i 
      = [ (u, this):rest 
          | this <- row (bs !! i) v (as !! i) lb 
          , let lb' = (replicate d 1 ++ map (+1) this) 
          , rest <- worker us vs ds lb' (i + 1)] 
    worker []     _      _      _  _ = [ [] ]
    worker (_:_)  []     _      _  _ = [ [] ]
    worker (_:_)  (_:_)  []     _  _ = [ [] ]
    -- weekly increasing lists of length @len@, pointwise at least @xs@, 
    -- maximum value @n@, minimum value @prev@.
    row :: Int -> Int -> Int -> [Int] -> [[Int]]
    row n len prev xxs = 
      if len == 0 
        then [[]] 
        else [ (j:js) | j <- [max x prev .. n], js <- row n (len-1) j xs ]
        where
          (x, xs) = fromJust (uncons xxs)

skewTableauWeight :: [(Int, [Int])] -> [Int]
skewTableauWeight skewT = [count i | i <- [1 .. m]]
  where
    (_, entries) = unzip skewT
    x = concat entries
    m = maximum x
    count i = sum [fromEnum (k == i) | k <- x]

isIncreasing :: [Int] -> Bool
isIncreasing s = 
  and (zipWith (<=) s (drop1 s))

-- _paths :: Int -> Seq Int -> Seq Int -> [[Seq Int]]
-- _paths n lambda mu = 
--   -- TODO: use same technique as macdonaldPolynomial: take only the compositions
--   -- which are partitions and use permuteMultiset
--   concatMap 
--     (skewGelfandTsetlinPatterns (DF.toList lambda) (DF.toList mu))
--       (compositions n (DF.sum lambda - DF.sum mu))

_paths :: Int -> Seq Int -> Seq Int -> [(Partition, [[(Seq Int, Seq Int)]])]
_paths n lambda mu =
  filter ((not . null) . snd) (map 
    (\nu -> let nu' = fromPartition nu 
                nu'' = nu' ++ replicate (n - length nu') 0
            in
      (
        nu''
      , map pairing (skewGelfandTsetlinPatterns lambda' (DF.toList mu) nu'')
      )
    ) 
    nus)
  where
    pairing lambdas = zip (drop1 lambdas) lambdas
    lambda' = DF.toList lambda
    nus = 
      filter ((<= n) . partitionWidth) $ 
        dominatedPartitions 
          (toPartitionUnsafe (lastSubPartition (DF.sum lambda - DF.sum mu) lambda'))

psi_lambda_mu :: forall a. (Eq a, AlgRing.C a) 
  => Seq Int -> Seq Int -> Spray a
psi_lambda_mu lambda mu = if S.null lambda
  then unitSpray
  else productOfSprays sprays
  where
    range = [1 .. lambda `S.index` 0]
    pair j = (
        1 + DF.sum (fmap (\k -> fromEnum (k == j)) lambda)
      , DF.sum (fmap (\k -> fromEnum (k == j)) mu)
      )
    pairs = filter (\(l, m) -> l == m) (map pair range)
    t = lone' 1
    sprays = map (\(_, m) -> AlgRing.one +> AlgAdd.negate (t m)) pairs

phi_lambda_mu :: forall a. (Eq a, AlgRing.C a) 
  => Seq Int -> Seq Int -> Spray a
phi_lambda_mu lambda mu = if S.null lambda
  then unitSpray
  else productOfSprays sprays
  where
    range = [1 .. lambda `S.index` 0]
    pair j = (
        DF.sum (fmap (\k -> fromEnum (k == j)) lambda)
      , 1 + DF.sum (fmap (\k -> fromEnum (k == j)) mu)
      )
    pairs = filter (\(l, m) -> l == m) (map pair range)
    t = lone' 1
    sprays = map (\(m, _) -> AlgRing.one +> AlgAdd.negate (t m)) pairs

_skewHallLittlewood :: forall a. (Eq a, AlgRing.C a) 
  => (Seq Int -> Seq Int -> Spray a) -> Int -> Seq Int -> Seq Int 
      -> SimpleParametricSpray a
_skewHallLittlewood f n lambda mu = 
  sumOfSprays (concatMap sprays paths)
  where
    paths = _paths n lambda mu
    allPairs = nub (concat (concat (snd (unzip paths))))
    psis = 
      HM.fromList 
        (map (\pair -> (pair, uncurry f pair)) allPairs)
    dropTrailingZeros = S.dropWhileR (== 0)
    sprays (nu, listsOfPairs) =
      let  
        sprays' = 
          [productOfSprays [psis HM.! pair | pair <- pairs] 
            | pairs <- listsOfPairs]
        listOfPowers = 
          [Powers expnts (S.length expnts) | 
            compo <- permuteMultiset nu, 
            let expnts = dropTrailingZeros (S.fromList compo)]
        in
        [
          HM.singleton powers spray
          | spray <- sprays', powers <- listOfPowers
        ]

-- skewHallLittlewoodP :: forall a. (Eq a, AlgRing.C a) 
--   => Int -> Seq Int -> Seq Int -> SimpleParametricSpray a
-- skewHallLittlewoodP n lambda mu = 
--   sumOfSprays [productOfSprays $ sprays path | path <- paths]
--   where
--     paths = _paths n lambda mu
--     lones = [lone' i | i <- [1 .. n]]
--     sprays nu = 
--       [psi_lambda_mu next_nu_i nu_i *^ lone_i (DF.sum next_nu_i - DF.sum nu_i)
--         | (next_nu_i, nu_i, lone_i) <- zip3 (drop 1 nu) nu lones]

-- skewHallLittlewoodQ :: forall a. (Eq a, AlgRing.C a) 
--   => Int -> Seq Int -> Seq Int -> SimpleParametricSpray a
-- skewHallLittlewoodQ n lambda mu = 
--   sumOfSprays [productOfSprays $ sprays path | path <- paths]
--   where
--     paths = _paths n lambda mu
--     lones = [lone' i | i <- [1 .. n]]
--     sprays nu = 
--       [phi_lambda_mu next_nu_i nu_i *^ lone_i (DF.sum next_nu_i - DF.sum nu_i)
--         | (next_nu_i, nu_i, lone_i) <- zip3 (drop1 nu) nu lones]

skewHallLittlewoodP :: forall a. (Eq a, AlgRing.C a) 
  => Int -> Seq Int -> Seq Int -> SimpleParametricSpray a
skewHallLittlewoodP = _skewHallLittlewood psi_lambda_mu

skewHallLittlewoodQ :: forall a. (Eq a, AlgRing.C a) 
  => Int -> Seq Int -> Seq Int -> SimpleParametricSpray a
skewHallLittlewoodQ = _skewHallLittlewood phi_lambda_mu

charge :: Seq Int -> Int
charge w = if l == 0 || n == 1 then 0 else DF.sum indices' + charge w'
  where
    l = S.length w
    n = DF.maximum w
    (positions', indices') = 
      go 1 (S.singleton (fromJust $ S.elemIndexL 1 w)) (S.singleton 0)
    w' = DF.foldr S.deleteAt w (S.sort positions')
    go :: Int -> Seq Int -> Seq Int -> (Seq Int, Seq Int)
    go r positions indices 
      | r == n    = (positions, indices)
      | otherwise = go (r+1) (positions |> pos') (indices |> index')
        where
          pos = positions `S.index` (r-1)
          index = indices `S.index` (r-1)
          v = S.drop (pos+1) w
          rindex = S.elemIndexL (r+1) v
          (pos', index') = 
            if isJust rindex
              then (1 + pos + fromJust rindex, index)
              else (fromJust (S.elemIndexL (r+1) w), index + 1)

_kostkaFoulkesPolynomial :: 
  (Eq a, AlgRing.C a) => Partition -> Partition -> Spray a
_kostkaFoulkesPolynomial lambda mu = 
  if sum lambda == sum mu 
    then sumOfSprays sprays
    else zeroSpray
  where
    tableaux = semiStandardTableauxWithGivenShapeAndWeight lambda mu
    mm = lone' 1 
    sprays =
      map (mm . charge . (mconcat . (map S.reverse))) tableaux 

b_lambda :: (Eq a, AlgRing.C a) => Partition -> Spray a
b_lambda lambda = productOfSprays sprays
  where
    table = [sum [fromEnum (k == j) | k <- lambda] | j <- nub lambda]
    sprays = map phi table
      where
        phi r = productOfSprays 
                [AlgRing.one +> AlgAdd.negate (lone' 1 i) | i <- [1 .. r]]

_transitionMatrixHallLittlewoodSchur :: 
  (Eq a, AlgRing.C a) => Char -> Int -> Map Partition (Map Partition (Spray a))
_transitionMatrixHallLittlewoodSchur which weight = 
  DM.fromDistinctDescList $ if which == 'P' 
    then zip lambdas [maps i | i <- rg]
    else zip lambdas 
              [DM.mapWithKey (\lambda c -> b_lambda lambda ^*^ c) (maps i) | i <- rg]
  where
    lambdas = reverse (map fromPartition (partitions weight))
    rg = [1 .. length lambdas]
    kfs = map f lambdas
    f kappa = 
      map (\mu -> _kostkaFoulkesPolynomial kappa mu)
          lambdas 
    matrix = inverseUnitTriangularMatrix (fromLists kfs)
    maps i = DM.filter (not . isZeroSpray) 
          (DM.fromDistinctDescList (zip lambdas (V.toList (getRow i matrix))))

_hallLittlewoodPolynomialsInSchurBasis :: 
  (Eq a, AlgRing.C a) => Char -> Partition -> Map Partition (Spray a)
_hallLittlewoodPolynomialsInSchurBasis which lambda = 
  if which == 'P'
    then coeffs
    else DM.map ((^*^) (b_lambda lambda)) coeffs
  where
    weight = sum lambda
    lambdas = 
      reverse $ filter (<= lambda) (map fromPartition (partitions weight))
    kfs = map f lambdas
    f kappa = 
      map (\mu -> _kostkaFoulkesPolynomial kappa mu) 
          lambdas -- (dominatedPartitions kappa)
    matrix = inverseUnitTriangularMatrix (fromLists kfs)
    coeffs = DM.filter (not . isZeroSpray) 
          (DM.fromDistinctDescList (zip lambdas (V.toList (getRow 1 matrix))))

_e :: AlgRing.C a => MCP.Partition -> a -> a
_e lambda alpha = 
  alpha * fromIntegral (_n (dualPartition lambda)) - fromIntegral (_n lambda)
  where
    _n mu = sum (zipWith (P.*) [0 .. ] (fromPartition mu))

_eSymbolic :: (Eq a, AlgRing.C a) => MCP.Partition -> Spray a 
_eSymbolic lambda = 
  _n (dualPartition lambda) .^ alpha <+ fromIntegral (- _n lambda)
  where
    alpha = lone 1
    _n mu = sum (zipWith (P.*) [0 .. ] (fromPartition mu))

_inverseKostkaMatrix :: 
  forall a. (Eq a, AlgField.C a) 
  => Int -> Int -> a -> Char -> (Matrix a, [Partition])
_inverseKostkaMatrix n weight alpha which = 
  (inverseTriangularMatrix (fromLists (map row lambdas)), lambdas)
  where
    kostkaNumbers = _kostkaNumbers n weight alpha which
    lambdas = reverse $ DM.keys kostkaNumbers
    msCombo lambda = kostkaNumbers DM.! lambda
    row lambda = 
      map (flip (DM.findWithDefault AlgAdd.zero) (msCombo lambda)) lambdas

_kostkaNumbers :: 
  forall a. (AlgField.C a) 
  => Int -> Int -> a -> Char -> Map Partition (Map Partition a)
_kostkaNumbers nv weight alpha which = kostkaMatrix'
  where
    coeffsP = DM.fromDistinctDescList 
      [(kappa, recip (jackCoeffP kappa alpha))| kappa <- lambdas']
    coeffsC = DM.fromDistinctDescList 
      [(kappa, jackCoeffC kappa alpha / jackCoeffP kappa alpha) 
        | kappa <- lambdas'] 
    coeffsQ = DM.fromDistinctDescList 
      [(kappa, jackCoeffQ kappa alpha / jackCoeffP kappa alpha) 
        | kappa <- lambdas']    
    kostkaMatrix = DM.mapKeys fromPartition (rec (length lambdas))
    kostkaMatrix' = case which of
      'J' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsP DM.! kappa)) m) 
                            kostkaMatrix
      'P' -> kostkaMatrix
      'C' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsC DM.! kappa)) m) 
                            kostkaMatrix
      'Q' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsQ DM.! kappa)) m) 
                            kostkaMatrix
      _   -> error "_kostkaNumbers: should not happen."
    mu_r_plus :: 
      Seq Int -> (Int, Int) -> Int -> (MCP.Partition, (Int, Int), Int)
    mu_r_plus mu pair@(i, j) r = 
      (
        MCP.Partition $ 
          DF.toList $ S.dropWhileR (== 0) $ S.reverse $ S.sort $ 
            S.adjust' ((P.+) r) i (S.adjust' (subtract r) j mu)
        , pair
        , r
      )
    lambdas = reverse $ 
      filter (\part -> partitionWidth part <= nv) (partitions weight)
    lambdas' = map fromPartition lambdas
    rec :: Int -> Map MCP.Partition (Map Partition a)
    rec n = if n == 1
      then DM.singleton (MCP.Partition [weight]) 
                        (DM.singleton [weight] AlgRing.one)
      else DM.insert mu (DM.singleton mu' AlgRing.one) 
            (
              DM.fromDistinctDescList 
              [(
                  kappa
                , DM.insert mu' (newColumn DM.! kappa) (previous DM.! kappa)
               ) | kappa <- kappas]
            ) 
      where
        previous = rec (n - 1)
        parts = take n lambdas
        (kappas, mu) = fromJust (unsnoc parts)
        _e_mu_alpha = _e mu alpha
        mu' = fromPartition mu
        mu'' = S.fromList mu'
        l = S.length mu''
        pairs = [(i, j) | i <- [0 .. l-2], j <- [i+1 .. l-1]]
        triplets = [mu_r_plus mu'' (i, j) r 
                    | (i, j) <- pairs, r <- [1 .. S.index mu'' j]]
        newColumn = 
          DM.fromDistinctDescList [(kappa, f kappa) | kappa <- kappas]
        f kappa = AlgAdd.sum xs 
          where
            previousRow = previous DM.! kappa
            triplets' = filter ((dominates kappa) . fst3) triplets
            ee = _e kappa alpha - _e_mu_alpha
            xs = [
              fromIntegral (S.index mu'' i P.- S.index mu'' j P.+ 2 P.* r) 
              * (previousRow DM.! (fromPartition nu)) / ee 
              | (nu, (i, j), r) <- triplets'
              ]

_symbolicKostkaNumbers :: 
  forall a. (Eq a, AlgField.C a) 
  => Int -> Int -> Char -> Map Partition (Map Partition (RatioOfSprays a))
_symbolicKostkaNumbers nv weight which = kostkaMatrix'
  where
    coeffsP = DM.fromDistinctDescList 
      [(kappa, asRatioOfSprays (jackSymbolicCoeffPinv kappa))
        | kappa <- lambdas']
    coeffsC = DM.fromDistinctDescList 
      [(
          kappa
        , (jackSymbolicCoeffPinv kappa :: Spray a) *> jackSymbolicCoeffC kappa
       ) | kappa <- lambdas']    
    coeffsQ = DM.fromDistinctDescList 
      [(
          kappa
        , jackSymbolicCoeffPinv kappa %//% jackSymbolicCoeffQinv kappa
       ) | kappa <- lambdas']    
    kostkaMatrix = DM.mapKeys fromPartition (rec (length lambdas))
    kostkaMatrix' = case which of
      'J' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsP DM.! kappa)) m) 
              kostkaMatrix
      'P' -> kostkaMatrix
      'C' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsC DM.! kappa)) m) 
              kostkaMatrix
      'Q' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsQ DM.! kappa)) m) 
              kostkaMatrix
      _   -> error "_symbolicKostkaNumbers: should not happen."
    mu_r_plus :: 
      Seq Int -> (Int, Int) -> Int -> (MCP.Partition, (Int, Int), Int)
    mu_r_plus mu pair@(i, j) r = 
      (
        MCP.Partition $ 
          DF.toList $ S.dropWhileR (== 0) $ S.reverse $ S.sort $ 
            S.adjust' ((P.+) r) i (S.adjust' (subtract r) j mu)
        , pair
        , r
      )
    lambdas = reverse $ 
      filter (\part -> partitionWidth part <= nv) (partitions weight)
    lambdas' = map fromPartition lambdas
    rec :: Int -> Map MCP.Partition (Map Partition (RatioOfSprays a))
    rec n = if n == 1
      then DM.singleton (MCP.Partition [weight]) 
                        (DM.singleton [weight] unitRatioOfSprays)
      else DM.insert mu (DM.singleton mu' unitRatioOfSprays) 
        (
          DM.fromDistinctDescList 
          [
            ( 
              kappa
            , DM.insert mu' (newColumn DM.! kappa) (previous DM.! kappa)
            ) 
            | kappa <- kappas
          ]
        ) 
      where
        previous = rec (n - 1)
        parts = take n lambdas
        (kappas, mu) = fromJust (unsnoc parts)
        _eSymbolic_mu = _eSymbolic mu
        mu' = fromPartition mu
        mu'' = S.fromList mu'
        l = S.length mu''
        pairs = [(i, j) | i <- [0 .. l-2], j <- [i+1 .. l-1]]
        triplets = [mu_r_plus mu'' (i, j) r 
                    | (i, j) <- pairs, r <- [1 .. S.index mu'' j]]
        newColumn = 
          DM.fromDistinctDescList [(kappa, f kappa) | kappa <- kappas]
        f kappa = AlgAdd.sum xs 
          where
            previousRow = previous DM.! kappa
            triplets' = filter ((dominates kappa) . fst3) triplets
            ee = _eSymbolic kappa - _eSymbolic_mu
            xs = [
              (
                (S.index mu'' i P.- S.index mu'' j P.+ 2 P.* r) 
                .^ (previousRow DM.! (fromPartition nu)) 
              ) %/% ee 
              | (nu, (i, j), r) <- triplets'
              ]

_inverseSymbolicKostkaMatrix :: 
  forall a. (Eq a, AlgField.C a) 
  => Int -> Int -> Char -> (Matrix (RatioOfSprays a), [Partition])
_inverseSymbolicKostkaMatrix n weight which = 
--  (inverseTriangularMatrix (fromLists (map (\lambda -> map (row lambda) lambdas) lambdas)), lambdas)
  (
    inverseTriangularMatrix (fromLists [map (row mu) lambdas | mu <- lambdas])
  , lambdas
  )
  where
    kostkaNumbers = _symbolicKostkaNumbers n weight which
    lambdas = reverse $ DM.keys kostkaNumbers
    msCombo lambda = kostkaNumbers DM.! lambda
    row = flip (DM.findWithDefault zeroRatioOfSprays) . msCombo
    -- row lambda = 
    --   map (flip (DM.findWithDefault zeroRatioOfSprays) (msCombo lambda)) lambdas

inverseTriangularMatrix :: (Eq a, AlgField.C a) => Matrix a -> Matrix a
inverseTriangularMatrix mat = 
  if d == 1 then fromLists [[recip (getElem 1 1 mat)]] else invmat
  where
    d = nrows mat
    invminor = inverseTriangularMatrix (minorMatrix d d mat)
    lastColumn = V.init (getCol d mat)
    vectors = [
        (
          V.drop (i-1) (getRow i invminor)
        , V.drop (i-1) lastColumn
        )
        | i <- [1 .. d-1]
      ] 
    lastEntry = recip (getElem d d mat)
    newColumn = colVector (V.fromList 
        [AlgAdd.negate (lastEntry * V.foldl1 (AlgAdd.+) (V.zipWith (*) u v)) 
          | (u, v) <- vectors]
      )
    newRow = rowVector (V.snoc (V.replicate (d - 1) AlgAdd.zero) lastEntry)
    invmat = (invminor <|> newColumn) <-> newRow

inverseUnitTriangularMatrix :: (Eq a, AlgRing.C a) => Matrix a -> Matrix a
inverseUnitTriangularMatrix mat = 
  if d == 1 then mat else invmat
  where
    d = nrows mat
    invminor = inverseUnitTriangularMatrix (minorMatrix d d mat)
    lastColumn = V.init (getCol d mat)
    vectors = [
        (
          V.drop (i-1) (getRow i invminor)
        , V.drop (i-1) lastColumn
        )
        | i <- [1 .. d-1]
      ] 
    newColumn = colVector (V.fromList 
        [AlgAdd.negate (V.foldl1 (AlgAdd.+) (V.zipWith (*) u v)) 
          | (u, v) <- vectors]
      )
    newRow = rowVector (V.snoc (V.replicate (d - 1) AlgAdd.zero) AlgRing.one)
    invmat = (invminor <|> newColumn) <-> newRow

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

hookLengths :: AlgRing.C a => Partition -> a -> ([a], [a])
hookLengths lambda alpha = (lower, upper)
  where
    (i, j) = _ij lambda
    (lambda', lambdaConj') = _convParts lambda
    upper = zipWith (fup lambdaConj' lambda') i j
      where
        fup x y ii jj =
          x!!(jj-1) - fromIntegral ii + 
            alpha * (y!!(ii-1) - fromIntegral (jj - 1))
    lower = zipWith (flow lambdaConj' lambda') i j
      where
        flow x y ii jj =
          x!!(jj-1) - (fromIntegral $ ii - 1) + 
            alpha * (y!!(ii-1) - fromIntegral jj)

_productHookLengths :: AlgRing.C a => Partition -> a -> a
_productHookLengths lambda alpha = product lower * product upper
  where
    (lower, upper) = hookLengths lambda alpha

jackCoeffC :: AlgField.C a => Partition -> a -> a
jackCoeffC lambda alpha = 
  alpha^k * fromInteger (product [2 .. k]) * recip jlambda
  where
    k = fromIntegral (sum lambda)
    jlambda = _productHookLengths lambda alpha

jackCoeffP :: AlgField.C a => Partition -> a -> a
jackCoeffP lambda alpha = one / product lower
  where
    (lower, _) = hookLengths lambda alpha

jackCoeffQ :: AlgField.C a => Partition -> a -> a
jackCoeffQ lambda alpha = one / product upper
  where
    (_, upper) = hookLengths lambda alpha

symbolicHookLengthsProducts :: forall a. (Eq a, AlgRing.C a) 
  => Partition -> (Spray a, Spray a)
symbolicHookLengthsProducts lambda = (product lower, product upper)
  where
    alpha = lone 1 :: Spray a
    (i, j) = _ij lambda
    (lambda', lambdaConj') = _convParts lambda
    upper = zipWith (fup lambdaConj' lambda') i j
      where
        fup x y ii jj =
          (x!!(jj-1) - fromIntegral ii) +> 
            ((y!!(ii-1) - fromIntegral (jj - 1)) *^ alpha)
    lower = zipWith (flow lambdaConj' lambda') i j
      where
        flow x y ii jj =
          (x!!(jj-1) - fromIntegral (ii - 1)) +> 
            ((y!!(ii-1) - fromIntegral jj) *^ alpha)

symbolicHookLengthsProduct :: (Eq a, AlgRing.C a) => Partition -> Spray a
symbolicHookLengthsProduct lambda = lower ^*^ upper
  where
    (lower, upper) = symbolicHookLengthsProducts lambda

jackSymbolicCoeffC :: 
  forall a. (Eq a, AlgField.C a) => Partition -> RatioOfSprays a
jackSymbolicCoeffC lambda = 
  ((fromIntegral factorialk) *^ alpha^**^k) %:% jlambda
  where
    alpha      = lone 1 :: Spray a
    k          = sum lambda
    factorialk = product [2 .. k]
    jlambda    = symbolicHookLengthsProduct lambda

jackSymbolicCoeffPinv :: (Eq a, AlgField.C a) => Partition -> Spray a
jackSymbolicCoeffPinv lambda = lower 
  where 
    (lower, _) = symbolicHookLengthsProducts lambda

jackSymbolicCoeffQinv :: (Eq a, AlgField.C a) => Partition -> Spray a 
jackSymbolicCoeffQinv lambda = upper 
  where 
    (_, upper) = symbolicHookLengthsProducts lambda

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

_betaRatioOfSprays :: forall a. (Eq a, AlgField.C a)
  => Partition -> Partition -> Int -> RatioOfSprays a
_betaRatioOfSprays kappa mu k = 
  ((x ^*^ num1 ^*^ num2 ^*^ num3) %:% (den1 ^*^ den2 ^*^ den3))
  where
    mukm1 = mu !! (k-1)
    x = lone 1 :: Spray a
    u = zipWith 
        (
        \s kap -> 
          (fromIntegral $ k - s + 1) +> ((fromIntegral $ kap - mukm1) *^ x)
        )
        [1 .. k] kappa 
    v = zipWith 
        (
        \s m -> (fromIntegral $ k - s) +> ((fromIntegral $ m - mukm1) *^ x)
        )
        [1 .. k-1] mu 
    w = zipWith 
        (
        \s m -> (fromIntegral $ m - k) +> ((fromIntegral $ mukm1 - s) *^ x)
        )
        [1 .. mukm1-1] (_dualPartition mu)
    num1 = product u
    den1 = product $ map (\p -> p ^+^ x ^-^ unitSpray) u
    num2 = product $ map (\p -> p ^+^ x) v
    den2 = product v
    num3 = product $ map (\p -> p ^+^ x) w
    den3 = product w

_fromInt :: (AlgRing.C a, Eq a) => Int -> a
_fromInt k = k .^ AlgRing.one

skewSchurLRCoefficients :: Partition -> Partition -> DM.Map Partition Int
skewSchurLRCoefficients lambda mu = 
  DM.mapKeys fromPartition (_lrRule lambda' mu')
  where
    lambda' = MCP.Partition lambda
    mu'     = MCP.Partition mu

isSkewPartition :: Partition -> Partition -> Bool
isSkewPartition lambda mu = 
  _isPartition lambda && _isPartition mu && all (>= 0) (zipWith (-) lambda mu)

sprayToMap :: Spray a -> Map [Int] a
sprayToMap spray = 
  DM.fromList (HM.toList $ HM.mapKeys (DF.toList . exponents) spray) 

comboToSpray :: (Eq a, AlgRing.C a) => Map Partition a -> Spray a
comboToSpray combo = sumOfSprays 
  [ let part' = S.fromList part in HM.singleton (Powers part' (S.length part')) c 
    | (part, c) <- DM.toList combo ]
