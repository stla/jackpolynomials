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
                                                             , foldl1'
                                                             , uncons
                                                             )
import           Data.List.Extra                             ( 
                                                               unsnoc
                                                             , snoc
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
                                                               Seq
                                                             , (|>) 
                                                             , (<|)
                                                             , (><)
                                                             , Seq ( (:<|) )
                                                             )
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
                                                             , SimpleParametricSpray
                                                             , zeroSpray
                                                             , isZeroSpray
                                                             , lone, lone', unitSpray
                                                             , sumOfSprays
                                                             , productOfSprays
                                                             , FunctionLike (..)
                                                             )
import           Math.Combinat.Partitions.Integer            (
                                                               fromPartition
                                                             , dualPartition
                                                             , partitions
                                                             , dominates
                                                             , partitionWidth
                                                             , toPartitionUnsafe
                                                             , dropTailingZeros
                                                             )
import qualified Math.Combinat.Partitions.Integer            as MCP
import           Math.Combinat.Partitions.Skew               (
                                                               SkewPartition
                                                             , mkSkewPartition
                                                             , skewPartitionElements
                                                             )
import           Math.Combinat.Tableaux.GelfandTsetlin       (
                                                                GT
                                                              , kostkaGelfandTsetlinPatterns
                                                             )
import           Math.Combinat.Tableaux.LittlewoodRichardson ( _lrRule )

type Partition = [Int]

gtPatternDiagonals :: GT -> (Int, [MCP.Partition])
gtPatternDiagonals pattern = (corner, [diagonal j | j <- [1 .. l]])
  where
    l = length pattern - 1
    corner = pattern !! l !! 0
    diagonal j = 
      (toPartitionUnsafe . dropTailingZeros) 
        [pattern !! r !! c | (r, c) <- zip [l-j .. l] [0 .. j]]

gtPatternToTableau :: GT -> [[Int]]
gtPatternToTableau pattern = DF.toList (go 0 startingTableau)
  where
    (corner, diagonals) = gtPatternDiagonals pattern
    diagonals' = toPartitionUnsafe [corner] : diagonals
    l = length diagonals - 1
    lambda = diagonals !! l
    m = partitionWidth lambda
    startingTableau = S.replicate m []
    zippedDiagonals = zip diagonals diagonals'
    skewPartition i = mkSkewPartition (zippedDiagonals !! i)
    go i tableau
      | i == 0 = go 1 (S.adjust' (++ replicate corner 1) 0 tableau)
      | i == l+2 = tableau
      | otherwise = 
          go (i+1) (growTableau (i+1) tableau (skewPartition (i-1)))
    growTableau :: Int -> Seq [Int] -> SkewPartition -> Seq [Int]
    growTableau j tableau skewPart =
      DF.foldr (\(i, _) -> S.adjust' (flip snoc j) (i-1)) tableau 
                (skewPartitionElements skewPart)

semiStandardTableauxWithGivenShapeAndWeight :: Partition -> Partition -> [[[Int]]]
semiStandardTableauxWithGivenShapeAndWeight lambda mu =
  map gtPatternToTableau (kostkaGelfandTsetlinPatterns lambda' mu')
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
  and [s !! i <= s !! (i+1) | i <- [0 .. length s - 2]]

isDecreasing :: Seq Int -> Bool
isDecreasing s = 
  and [s `S.index` i >= s `S.index` (i+1) | i <- [0 .. S.length s - 2]]

cartesianProduct :: Seq Int -> [Seq Int]
cartesianProduct (S.Empty) = []
cartesianProduct (i:<|is)
  | S.null is = [S.singleton j | j <- [i, i-1 .. 0]]
  | otherwise = [j <| s | j <- [i, i-1 .. 0], s <- previous]
    where
      previous = cartesianProduct is

horizontalStrip :: Seq Int -> Seq Int -> Bool
horizontalStrip lambda mu = all (`elem` [0, 1]) theta'
  where
    lambda' = S.fromList $ _dualPartition (DF.toList lambda)
    mu' = S.fromList $ _dualPartition (DF.toList mu)
    mu'' = mu' >< (S.replicate (S.length lambda' - S.length mu') 0)
    theta' = S.zipWith (-) lambda' mu''

columnStrictTableau :: [Seq Int] -> Bool
columnStrictTableau tableau = 
  and (zipWith horizontalStrip tableau tail_tableau)
  where tail_tableau = drop 1 tableau

_paths :: Int -> Seq Int -> Seq Int -> [[Seq Int]]
_paths n lambda mu = filter columnStrictTableau tableaux
  where
    mu' = mu >< (S.replicate (S.length lambda - S.length mu) 0)
    diffs = S.zipWith (-) lambda mu'
    grid = cartesianProduct diffs
    kappas = filter isDecreasing [S.zipWith (+) kappa mu' | kappa <- grid]
    combos = combinations 0 (length kappas - 1) (n-1)
      where
        combinations :: Int -> Int -> Int -> [[Int]]
        combinations a b m 
          | m == 0 = [[]]
          | m == 1 = [[i] | i <- [a .. b]]
          | otherwise = 
              [i : combo | i <- [a .. b], combo <- combinations i b (m-1)]
    tableaux = 
      map (\combo -> lambda : (map ((!!) kappas) combo) ++ [mu']) combos

psi_lambda_mu :: forall a. (Eq a, AlgRing.C a) 
  => Seq Int -> Seq Int -> Spray a
psi_lambda_mu lambda mu = productOfSprays sprays
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
phi_lambda_mu lambda mu = productOfSprays sprays
  where
    range = [1 .. lambda `S.index` 0]
    pair j = (
        DF.sum (fmap (\k -> fromEnum (k == j)) lambda)
      , 1 + DF.sum (fmap (\k -> fromEnum (k == j)) mu)
      )
    pairs = filter (\(l, m) -> l == m) (map pair range)
    t = lone' 1
    sprays = map (\(m, _) -> AlgRing.one +> AlgAdd.negate (t m)) pairs

skewHallLittlewoodP :: forall a. (Eq a, AlgRing.C a) 
  => Int -> Seq Int -> Seq Int -> SimpleParametricSpray a
skewHallLittlewoodP n lambda mu = 
  sumOfSprays [productOfSprays $ sprays (reverse path) | path <- paths]
  where
    paths = _paths n lambda mu
    lones = [lone' i | i <- [1 .. n]]
    sprays nu = 
      [psi_lambda_mu next_nu_i nu_i *^ lone_i (DF.sum next_nu_i - DF.sum nu_i)
        | (next_nu_i, nu_i, lone_i) <- zip3 (drop 1 nu) nu lones]

skewHallLittlewoodQ :: forall a. (Eq a, AlgRing.C a) 
  => Int -> Seq Int -> Seq Int -> SimpleParametricSpray a
skewHallLittlewoodQ n lambda mu = 
  sumOfSprays [productOfSprays $ sprays (reverse path) | path <- paths]
  where
    paths = _paths n lambda mu
    lones = [lone' i | i <- [1 .. n]]
    sprays nu = 
      [phi_lambda_mu next_nu_i nu_i *^ lone_i (DF.sum next_nu_i - DF.sum nu_i)
        | (next_nu_i, nu_i, lone_i) <- zip3 (drop 1 nu) nu lones]

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

isDominated :: Seq Int -> Seq Int -> Bool
isDominated mu lambda = 
  (MCP.Partition (DF.toList lambda)) `dominates` (MCP.Partition (DF.toList mu))

-- assumes sum lambda == sum mu 
ssytxWithGivenShapeAndContent :: Seq Int -> Seq Int -> [Seq (Seq Int)]
ssytxWithGivenShapeAndContent lambda mu = 
  if all (== 1) lambda 
    then if all (== 1) mu
      then [S.fromList [S.singleton i | i <- [1 .. S.length lambda]]]
      else []
    else if isDominated mu lambda
      then nub all_ssytx
      else []
  where
    dropTrailingZeros = S.dropWhileR (== 0)
    l = S.length lambda
    m = S.length mu
    mu' = dropTrailingZeros $ S.adjust' (subtract 1) (m-1) mu
    zippedKappas = 
      zip [0 ..] [S.adjust' (subtract 1) i lambda | i <- [0 .. l - 1]]
    all_ssytx = concatMap f zippedKappas
      where
        f (i, kappa) = 
           if isDecreasing kappa 
            then nub $ 
              map g (ssytxWithGivenShapeAndContent kappa' mu')
            else []
          where 
            kappa' = dropTrailingZeros kappa
            g ssyt = if i < S.length ssyt 
              then S.adjust' (|> m) i ssyt 
              else ssyt |> (S.singleton m)
            -- g ssyt = if i < length ssyt 
            --   then (element i .~ ssyt !! i |> l) ssyt 
            --   else ssyt ++ [S.singleton l]

_kostkaFoulkesPolynomial :: 
  (Eq a, AlgRing.C a) => Seq Int -> Seq Int -> Spray a
_kostkaFoulkesPolynomial lambda mu = 
  if DF.sum lambda == DF.sum mu 
    then sumOfSprays sprays
    else zeroSpray
  where
    tableaux = ssytxWithGivenShapeAndContent lambda mu
    mm = lone' 1 -- TODO: fix lone' 1 0 (= fromList [(Powers {exponents = fromList [0], nvariables = 1},1 % 1)])
    sprays = 
      map (mm . charge . ((foldl1' (S.><)) . (map S.reverse) . DF.toList)) 
            tableaux

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
      map (\mu -> _kostkaFoulkesPolynomial (S.fromList kappa) (S.fromList mu))
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
      map (\mu -> _kostkaFoulkesPolynomial (S.fromList kappa) (S.fromList mu)) 
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
