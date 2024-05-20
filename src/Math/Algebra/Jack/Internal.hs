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
  , inverseTriangularMatrix
  , _kostkaNumbers
  , _inverseKostkaMatrix
  , _symbolicKostkaNumbers
  , _inverseSymbolicKostkaMatrix
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
import           Data.List.Extra                             ( unsnoc )
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
import           Data.Maybe                                  ( fromJust )
import qualified Data.Sequence                               as S
import           Data.Tuple.Extra                            ( fst3 )
import qualified Data.Vector                                 as V
import           Math.Algebra.Hspray                         ( 
                                                               RatioOfSprays, (%:%), (%//%), (%/%)
                                                             , unitRatioOfSprays
                                                             , zeroRatioOfSprays
                                                             , asRatioOfSprays
                                                             , Spray, (.^)
                                                             , Powers (..)
                                                             , lone, unitSpray
                                                             , sumOfSprays
                                                             , FunctionLike (..)
                                                             )
import           Math.Combinat.Partitions.Integer            (
                                                               fromPartition
                                                             , dualPartition
                                                             , partitions
                                                             , dominates
                                                             , mkPartition
                                                             , partitionWidth
                                                             )
import qualified Math.Combinat.Partitions.Integer            as MCP
import           Math.Combinat.Tableaux.LittlewoodRichardson ( _lrRule )

type Partition = [Int]


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

_inverseKostkaMatrix :: forall a. (Eq a, AlgField.C a) => Int -> Int -> a -> Char -> (Matrix a, [Partition])
_inverseKostkaMatrix n weight alpha which = 
  (inverseTriangularMatrix (fromLists (map row lambdas)), lambdas)
  where
    kostkaNumbers = _kostkaNumbers n weight alpha which
    lambdas = reverse $ DM.keys kostkaNumbers
    msCombo lambda = kostkaNumbers DM.! lambda
    row lambda = map (flip (DM.findWithDefault AlgAdd.zero) (msCombo lambda)) lambdas

_kostkaNumbers :: forall a. (AlgField.C a) => Int -> Int -> a -> Char -> Map Partition (Map Partition a)
_kostkaNumbers nv weight alpha which = kostkaMatrix'
  where
    coeffsP = DM.fromList 
      [(kappa, recip (jackCoeffP kappa alpha))| kappa <- lambdas']
    coeffsC = DM.fromList 
      [(kappa, jackCoeffC kappa alpha / jackCoeffP kappa alpha) | kappa <- lambdas']    
    coeffsQ = DM.fromList 
      [(kappa, jackCoeffQ kappa alpha / jackCoeffP kappa alpha) | kappa <- lambdas']    
    kostkaMatrix = DM.mapKeys fromPartition (rec (length lambdas))
    kostkaMatrix' = case which of
      'J' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsP DM.! kappa)) m) kostkaMatrix
      'P' -> kostkaMatrix
      'C' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsC DM.! kappa)) m) kostkaMatrix
      'Q' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsQ DM.! kappa)) m) kostkaMatrix
      _   -> error "_kostkaNumbers: should not happen."
    mu_r_plus :: Partition -> (Int, Int) -> Int -> (MCP.Partition, (Int, Int), Int)
    mu_r_plus mu pair@(i, j) r = 
      (mkPartition (DF.toList $ S.reverse $ S.sort $ (S.adjust' ((P.+) r) i (S.adjust' (subtract r) j mu'))), pair, r)
      where
        mu' = S.fromList mu 
    lambdas = reverse $ filter (\part -> partitionWidth part <= nv) (partitions weight)
    lambdas' = map fromPartition lambdas
    rec :: Int -> Map MCP.Partition (Map Partition a)
    rec n = if n == 1
      then DM.singleton (MCP.Partition [weight]) (DM.singleton [weight] AlgRing.one)
      else DM.insert mu (DM.singleton mu' AlgRing.one) (DM.fromList [(kappa, DM.insert mu' (newColumn DM.! kappa) (previous DM.! kappa)) | kappa <- kappas]) --  DM.union previous (DM.singleton mu newColumn
      where
        previous = rec (n - 1)
        parts = take (n) lambdas
        (kappas, mu) = fromJust (unsnoc parts)
        mu' = fromPartition mu
        l = length mu'
        pairs = [(i, j) | i <- [0 .. l-2], j <- [i+1 .. l-1]]
        newColumn = DM.fromList [(kappa, f kappa) | kappa <- kappas]
        f kappa = AlgAdd.sum xs 
          where
            previousRow = previous DM.! kappa
            nus = 
              filter ((dominates kappa) . fst3) [mu_r_plus mu' (i, j) r | (i, j) <- pairs, r <- [1 .. mu' !! j]]
            ee = _e kappa alpha - _e mu alpha
            xs = [fromIntegral (mu' !! i P.- mu' !! j P.+ 2 P.* r) * (previousRow DM.! (fromPartition nu)) / ee | (nu, (i, j), r) <- nus]

_symbolicKostkaNumbers :: forall a. (Eq a, AlgField.C a) => Int -> Int -> Char -> Map Partition (Map Partition (RatioOfSprays a))
_symbolicKostkaNumbers nv weight which = kostkaMatrix'
  where
    coeffsP = DM.fromList 
      [(kappa, asRatioOfSprays (jackSymbolicCoeffPinv kappa))| kappa <- lambdas']
    coeffsC = DM.fromList 
      [(kappa, (jackSymbolicCoeffPinv kappa :: Spray a) *> jackSymbolicCoeffC kappa) | kappa <- lambdas']    
    coeffsQ = DM.fromList 
      [(kappa, jackSymbolicCoeffPinv kappa %//% jackSymbolicCoeffQinv kappa) | kappa <- lambdas']    
    kostkaMatrix = DM.mapKeys fromPartition (rec (length lambdas))
    kostkaMatrix' = case which of
      'J' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsP DM.! kappa)) m) kostkaMatrix
      'P' -> kostkaMatrix
      'C' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsC DM.! kappa)) m) kostkaMatrix
      'Q' -> DM.mapWithKey (\kappa m -> DM.map ((*) (coeffsQ DM.! kappa)) m) kostkaMatrix
      _   -> error "_kostkaNumbers: should not happen."
    mu_r_plus :: Partition -> (Int, Int) -> Int -> (MCP.Partition, (Int, Int), Int)
    mu_r_plus mu pair@(i, j) r = 
      (mkPartition (DF.toList $ S.reverse $ S.sort $ (S.adjust' ((P.+) r) i (S.adjust' (subtract r) j mu'))), pair, r)
      where
        mu' = S.fromList mu 
    lambdas = reverse $ filter (\part -> partitionWidth part <= nv) (partitions weight)
    lambdas' = map fromPartition lambdas
    rec :: Int -> Map MCP.Partition (Map Partition (RatioOfSprays a))
    rec n = if n == 1
      then DM.singleton (MCP.Partition [weight]) (DM.singleton [weight] unitRatioOfSprays)
      else DM.insert mu (DM.singleton mu' unitRatioOfSprays) (DM.fromList [(kappa, DM.insert mu' (newColumn DM.! kappa) (previous DM.! kappa)) | kappa <- kappas]) 
      where
        previous = rec (n - 1)
        parts = take (n) lambdas
        (kappas, mu) = fromJust (unsnoc parts)
        mu' = fromPartition mu
        l = length mu'
        pairs = [(i, j) | i <- [0 .. l-2], j <- [i+1 .. l-1]]
        newColumn = DM.fromList [(kappa, f kappa) | kappa <- kappas]
        f kappa = AlgAdd.sum xs 
          where
            previousRow = previous DM.! kappa
            nus = 
              filter ((dominates kappa) . fst3) [mu_r_plus mu' (i, j) r | (i, j) <- pairs, r <- [1 .. (mu' !! j P.+ mu' !! j) `div` 2]]
            ee = _eSymbolic kappa - _eSymbolic mu
            xs = [((mu' !! i P.- mu' !! j P.+ 2 P.* r) .^ (previousRow DM.! (fromPartition nu))) %/% ee | (nu, (i, j), r) <- nus]

_inverseSymbolicKostkaMatrix :: 
  forall a. (Eq a, AlgField.C a) => Int -> Int -> Char -> (Matrix (RatioOfSprays a), [Partition])
_inverseSymbolicKostkaMatrix n weight which = 
  (inverseTriangularMatrix (fromLists (map row lambdas)), lambdas)
  where
    kostkaNumbers = _symbolicKostkaNumbers n weight which
    lambdas = reverse $ DM.keys kostkaNumbers
    msCombo lambda = kostkaNumbers DM.! lambda
    row lambda = map (flip (DM.findWithDefault zeroRatioOfSprays) (msCombo lambda)) lambdas










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
