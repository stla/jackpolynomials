{-|
Module      : Math.Algebra.Jack.SymmetricPolynomials
Description : Some utilities for Jack polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

A Jack polynomial can have a very long expression in the canonical basis. 
A considerably shorter expression is obtained by writing the polynomial as 
a linear combination of the monomial symmetric polynomials instead, which is 
always possible since Jack polynomials are symmetric. This is the initial 
motivation of this module. But now it contains more stuff dealing with 
symmetric polynomials.
-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Math.Algebra.SymmetricPolynomials
  ( 
  -- * Checking symmetry
    isSymmetricSpray
  -- * Classical symmetric polynomials
  , msPolynomial
  , psPolynomial
  , cshPolynomial
  -- * Decomposition of symmetric polynomials
  , msCombination
  , psCombination
  , psCombination'
  -- * Printing symmetric polynomials
  , prettySymmetricNumSpray
  , prettySymmetricQSpray
  , prettySymmetricQSpray'
  , prettySymmetricParametricQSpray
  -- * Operators on the space of symmetric polynomials
  , laplaceBeltrami
  , calogeroSutherland
  -- * Hall inner product of symmetric polynomials
  , hallInnerProduct
  , hallInnerProduct'
  , hallInnerProduct''
  , hallInnerProduct'''
  , hallInnerProduct''''
  , symbolicHallInnerProduct
  , symbolicHallInnerProduct'
  , symbolicHallInnerProduct''
  , test
  ) where
import           Prelude hiding ( fromIntegral, fromRational )
import qualified Algebra.Additive                 as AlgAdd
import           Algebra.Field                    ( fromRational )
import qualified Algebra.Field                    as AlgField
import qualified Algebra.Module                   as AlgMod
import qualified Algebra.Ring                     as AlgRing
import           Algebra.ToInteger                ( fromIntegral )
import qualified Data.Foldable                    as DF
import qualified Data.HashMap.Strict              as HM
import           Data.List                        ( foldl1', nub )
import           Data.List.Extra                  ( unsnoc )
import           Data.Map.Strict                  ( 
                                                    Map
                                                  , unionsWith
                                                  , insert
                                                  )
import qualified Data.Map.Strict                  as DM
import           Data.Maybe                       ( fromJust )
import           Data.Map.Merge.Strict            ( 
                                                    merge
                                                  , dropMissing
                                                  , zipWithMatched 
                                                  )
import           Data.Ratio                       ( (%) )
import           Data.Sequence                    ( 
                                                    Seq
                                                  , (|>) 
                                                  )
import qualified Data.Sequence                    as S
import           Data.Tuple.Extra                 ( second )
import           Math.Algebra.Hspray              (
                                                    FunctionLike (..)
                                                  , (/^)
                                                  , Spray
                                                  , Powers (..)
                                                  , QSpray
                                                  , QSpray'
                                                  , ParametricQSpray
                                                  , lone
                                                  , qlone
                                                  , lone'
                                                  , fromList
                                                  , getCoefficient
                                                  , getConstantTerm
                                                  , isConstant
                                                  , (%//%)
                                                  , RatioOfSprays (..)
                                                  , prettyRatioOfQSpraysXYZ
                                                  , showNumSpray
                                                  , showQSpray
                                                  , showQSpray'
                                                  , showSpray
                                                  , toList
                                                  , zeroSpray
                                                  , unitSpray
                                                  , productOfSprays
                                                  , sumOfSprays
                                                  , constantSpray
                                                  )
import           Math.Algebra.Jack.Internal       ( Partition , _isPartition )
import           Math.Combinat.Compositions       ( compositions1 )
import           Math.Combinat.Partitions.Integer ( 
                                                    fromPartition
                                                  , mkPartition
                                                  , partitions 
                                                  )
import           Math.Combinat.Permutations       ( permuteMultiset )

-- | Monomial symmetric polynomial
--
-- >>> putStrLn $ prettySpray' (msPolynomial 3 [2, 1])
-- (1) x1^2.x2 + (1) x1^2.x3 + (1) x1.x2^2 + (1) x1.x3^2 + (1) x2^2.x3 + (1) x2.x3^2
msPolynomial :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
msPolynomial n lambda
  | n < 0                     = 
      error "msPolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "msPolynomial: invalid partition."
  | llambda > n               = zeroSpray
  | otherwise                 = fromList $ zip permutations coefficients
    where
      llambda      = length lambda
      permutations = permuteMultiset (lambda ++ replicate (n-llambda) 0)
      coefficients = repeat AlgRing.one

-- | Checks whether a spray defines a symmetric polynomial; this is useless for 
-- Jack polynomials because they always are symmetric, but this module contains 
-- everything needed to build this function which can be useful in another context
isSymmetricSpray :: (AlgRing.C a, Eq a) => Spray a -> Bool
isSymmetricSpray spray = spray == spray' 
  where
    assocs = msCombination' spray
    n      = numberOfVariables spray
    spray' = foldl1' (^+^) 
      (
        map (\(lambda, x) -> x *^ msPolynomial n lambda) assocs
      )

-- | Symmetric polynomial as a linear combination of monomial symmetric polynomials
msCombination :: AlgRing.C a => Spray a -> Map Partition a
msCombination spray = DM.fromList (msCombination' spray)

msCombination' :: AlgRing.C a => Spray a -> [(Partition, a)]
msCombination' spray = 
  map (\lambda -> (lambda, getCoefficient lambda spray)) lambdas
  where
    lambdas = nub $ map (fromPartition . mkPartition . fst) (toList spray)

-- helper function for the showing stuff
makeMSpray :: (Eq a, AlgRing.C a) => Spray a -> Spray a
makeMSpray = fromList . msCombination'

-- show symmetric monomial like M[3,2,1]
showSymmetricMonomials :: [Seq Int] -> [String]
showSymmetricMonomials = map showSymmetricMonomial
  where
    showSymmetricMonomial :: Seq Int -> String
    showSymmetricMonomial lambda = 'M' : show (DF.toList lambda)

-- | Prints a symmetric spray as a linear combination of monomial symmetric polynomials
--
-- >>> putStrLn $ prettySymmetricNumSpray $ schurPol' 3 [3, 1, 1]
-- M[3,1,1] + M[2,2,1]
prettySymmetricNumSpray :: 
  (Num a, Ord a, Show a, AlgRing.C a) => Spray a -> String
prettySymmetricNumSpray spray = 
  showNumSpray showSymmetricMonomials show mspray
  where
    mspray = makeMSpray spray

-- | Prints a symmetric spray as a linear combination of monomial symmetric polynomials
--
-- >>> putStrLn $ prettySymmetricQSpray $ jackPol' 3 [3, 1, 1] 2 'J'
-- 42*M[3,1,1] + 28*M[2,2,1]
prettySymmetricQSpray :: QSpray -> String
prettySymmetricQSpray spray = showQSpray showSymmetricMonomials mspray
  where
    mspray = makeMSpray spray

-- | Same as `prettySymmetricQSpray` but for a `QSpray'` symmetric spray
prettySymmetricQSpray' :: QSpray' -> String
prettySymmetricQSpray' spray = showQSpray' showSymmetricMonomials mspray
  where
    mspray = makeMSpray spray

-- | Prints a symmetric parametric spray as a linear combination of monomial 
-- symmetric polynomials
--
-- >>> putStrLn $ prettySymmetricParametricQSpray ["a"] $ jackSymbolicPol' 3 [3, 1, 1] 'J'
-- { [ 4*a^2 + 10*a + 6 ] }*M[3,1,1] + { [ 8*a + 12 ] }*M[2,2,1]
prettySymmetricParametricQSpray :: [String] -> ParametricQSpray -> String
prettySymmetricParametricQSpray letters spray = 
  showSpray (prettyRatioOfQSpraysXYZ letters) ("{ ", " }") 
            showSymmetricMonomials mspray
  where
    mspray = makeMSpray spray

-- | Laplace-Beltrami operator on the space of homogeneous symmetric polynomials;
-- neither symmetry and homogeneity are checked
laplaceBeltrami :: (Eq a, AlgField.C a) => a -> Spray a -> Spray a
laplaceBeltrami alpha spray = 
  if isConstant spray 
    then zeroSpray 
    else alpha' *^ spray1 ^+^ spray2
  where
    alpha' = alpha AlgField./ AlgRing.fromInteger 2
    n = numberOfVariables spray
    range = [1 .. n]
    dsprays = map (`derivative` spray) range
    op1 i = lone' i 2 ^*^ derivative i (dsprays !! (i-1))
    spray1 = AlgAdd.sum (map op1 range)
    spray2 = _numerator $ AlgAdd.sum 
              [(lone' i 2 ^*^ dsprays !! (i-1)) %//% (lone i ^-^ lone j)
                | i <- range, j <- range, i /= j]

-- | Calogero-Sutherland operator on the space of homogeneous symmetric polynomials;
-- neither symmetry and homogeneity are checked
calogeroSutherland :: (Eq a, AlgField.C a) => a -> Spray a -> Spray a
calogeroSutherland alpha spray = 
  if isConstant spray 
    then zeroSpray
    else halfSpray $ alpha *^ spray1 ^+^ spray2
  where
    halfSpray p = p /^ AlgRing.fromInteger 2
    n = numberOfVariables spray
    range = [1 .. n]
    dsprays = map (`derivative` spray) range
    op0 p i = lone i ^*^ derivative i p 
    op1 p i = op0 (op0 p i) i
    spray1 = AlgAdd.sum (map (op1 spray) range)
    spray2 = _numerator $ AlgAdd.sum 
      [let (xi, xj, dxi, dxj) = 
            (lone i, lone j, dsprays !! (i-1), dsprays !! (j-1)) in 
          (xi ^+^ xj) ^*^ (xi ^*^ dxi ^-^ xj ^*^ dxj) %//% (xi ^-^ xj)
       | i <- range, j <- [i+1 .. n]]

-- | Power sum polynomial
--
-- >>> putStrLn $ prettyQSpray (psPolynomial 3 [2, 1])
-- x^3 + x^2.y + x^2.z + x.y^2 + x.z^2 + y^3 + y^2.z + y.z^2 + z^3
psPolynomial :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
psPolynomial n lambda
  | n < 0                     = 
      error "psPolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "psPolynomial: invalid partition."
  | null lambda               = unitSpray
  | llambda > n               = zeroSpray
  | otherwise                 = productOfSprays sprays
    where
      llambda = length lambda
      sprays = [HM.fromList $ [f i k | i <- [1 .. n]] | k <- lambda]
      f j k = (Powers expts j, AlgRing.one)
        where
          expts = S.replicate (j-1) 0 |> k

eLambdaMu :: Partition -> Partition -> Rational
eLambdaMu lambda mu 
  | ellLambda < ellMu = 0
  | otherwise = if even (ellLambda - ellMu) 
      then sum xs 
      else - sum xs
  where
    ellLambda = length lambda
    ellMu     = length mu
    compos = compositions1 ellMu ellLambda
    lambdaPerms = permuteMultiset lambda
    sequencesOfPartitions = filter (not . null)
      [partitionSequences perm compo 
        | perm <- lambdaPerms, compo <- compos]
    xs = [eMuNus nus | nus <- sequencesOfPartitions]
    ----
    partitionSequences :: [Int] -> [Int] -> [Partition]
    partitionSequences kappa compo = if test then nus else []
      where
        headOfCompo = fst $ fromJust (unsnoc compo)
        starts = scanl (+) 0 headOfCompo 
        ends   = zipWith (+) starts compo
        nus = [ 
                [ kappa !! k | k <- [starts !! i .. ends !! i - 1] ] 
                | i <- [0 .. length compo - 1]
              ]
        nuWeights = [sum nu | nu <- nus]
        decreasing ys = 
          and [ys !! i >= ys !! (i+1) | i <- [0 .. length ys - 2]]
        test = and (zipWith (==) mu nuWeights) && all decreasing nus
    ---- 
    eMuNus :: [Partition] -> Rational
    eMuNus nus = product toMultiply
      where
        w :: Int -> Partition -> Rational
        w k nu = 
          let table = [sum [fromEnum (i == j) | i <- nu] | j <- nub nu] in
          (toInteger $ k * factorial (length nu - 1)) % 
            (toInteger $ product (map factorial table))
        factorial n = product [2 .. n]
        toMultiply = zipWith w mu nus

-- | monomial symmetric polynomial as a linear combination of 
-- power sum polynomials
mspInPSbasis :: Partition -> Map Partition Rational
mspInPSbasis kappa = DM.fromList (zipWith f weights lambdas)
  where
    parts = map fromPartition (partitions (sum kappa))
    (weights, lambdas) = unzip $ filter ((/= 0) . fst) 
      [(eLambdaMu kappa lambda, lambda) | lambda <- parts]
    f weight lambda = 
      (lambda, weight / toRational (zlambda lambda))

-- | the factor in the Hall inner product
zlambda :: Partition -> Int
zlambda lambda = p
  where
    parts = nub lambda
    table = [sum [fromEnum (k == j) | k <- lambda] | j <- parts]
    p =  
      product [factorial mj * part^mj | (part, mj) <- zip parts table]
    factorial n = product [2 .. n]

-- | symmetric polynomial as a linear combination of power sum polynomials
_psCombination :: 
  forall a. (Eq a, AlgRing.C a) => (a -> Rational -> a) -> Spray a -> Map Partition a
_psCombination func spray =
  if constantTerm == AlgAdd.zero 
    then psMap
    else insert [] constantTerm psMap
  where
    constantTerm = getConstantTerm spray
    assocs = msCombination' (spray <+ (AlgAdd.negate constantTerm))
    f :: (Partition, a) -> [(Partition, a)] 
    f (lambda, coeff) = 
      map (second (func coeff)) (DM.toList psCombo)
      where
        psCombo = mspInPSbasis lambda :: Map Partition Rational
    psMap = DM.filter (/= AlgAdd.zero) 
            (unionsWith (AlgAdd.+) (map (DM.fromList . f) assocs))

-- | Symmetric polynomial as a linear combination of power sum polynomials. 
-- Symmetry is not checked.
psCombination :: 
  forall a. (Eq a, AlgField.C a) => Spray a -> Map Partition a
psCombination = 
  _psCombination (\coef r -> coef AlgRing.* fromRational r)

-- | Symmetric polynomial as a linear combination of power sum polynomials. 
-- Same as @psCombination@ but with other constraints on the base ring of the spray.
psCombination' :: 
  forall a. (Eq a, AlgMod.C Rational a, AlgRing.C a) 
  => Spray a -> Map Partition a
psCombination' = _psCombination (flip (AlgMod.*>))

-- | the Hall inner product with parameter
_hallInnerProduct :: 
  forall a b. (AlgRing.C b, AlgRing.C a)
  => (Spray b -> Map Partition b)
  -> (a -> b -> b)
  -> Spray b   -- ^ spray
  -> Spray b   -- ^ spray
  -> a         -- ^ parameter
  -> b 
_hallInnerProduct psCombinationFunc multabFunc spray1 spray2 alpha = 
  AlgAdd.sum $ DM.elems
    (merge dropMissing dropMissing (zipWithMatched f) psCombo1 psCombo2)
  where
    psCombo1 = psCombinationFunc spray1 :: Map Partition b
    psCombo2 = psCombinationFunc spray2 :: Map Partition b
    zlambda' :: Partition -> a
    zlambda' lambda = fromIntegral (zlambda lambda) 
      AlgRing.* alpha AlgRing.^ (toInteger $ length lambda)
    f :: Partition -> b -> b -> b
    f lambda coeff1 coeff2 = 
      multabFunc (zlambda' lambda) (coeff1 AlgRing.* coeff2)

-- | Hall inner product with parameter. It makes sense only for symmetric sprays,
-- and the symmetry is not checked. 
hallInnerProduct :: 
  forall a. (Eq a, AlgField.C a)
  => Spray a   -- ^ spray
  -> Spray a   -- ^ spray
  -> a         -- ^ parameter
  -> a 
hallInnerProduct = _hallInnerProduct psCombination (AlgRing.*)

-- | Hall inner product with parameter. Same as @hallInnerProduct@ but 
-- with other constraints on the base ring of the sprays.
hallInnerProduct' :: 
  forall a. (Eq a, AlgMod.C Rational a, AlgRing.C a)
  => Spray a   -- ^ spray
  -> Spray a   -- ^ spray
  -> a         -- ^ parameter
  -> a 
hallInnerProduct' = _hallInnerProduct psCombination' (AlgRing.*)

-- | Hall inner product with parameter. Same as @hallInnerProduct@ but 
-- with other constraints on the base ring of the sprays. It is applicable 
-- to @Spray Int@ sprays.
hallInnerProduct'' :: 
  forall a. (Real a)
  => Spray a   -- ^ spray
  -> Spray a   -- ^ spray
  -> a         -- ^ parameter
  -> Rational 
hallInnerProduct'' spray1 spray2 alpha = 
  _hallInnerProduct 
    (_psCombination (*)) (*) qspray1 qspray2 (toRational alpha)
  where
    asQSpray :: Spray a -> QSpray
    asQSpray = HM.map toRational
    qspray1 = asQSpray spray1
    qspray2 = asQSpray spray2

-- | Hall inner product with parameter for parametric sprays, because the
-- type of the parameter in @hallInnerProduct@ is strange. For example, a
-- @ParametricQSpray@ spray is a @Spray RatioOfQSprays@ spray, and it makes
-- more sense to compute the Hall product with a @Rational@ parameter then 
-- to compute the Hall product with a @RatioOfQSprays@ parameter.
--
-- >>> import Math.Algebra.Jack.SymmetricPolynomials
-- >>> import Math.Algebra.JackSymbolicPol
-- >>> import Math.Algebra.Hspray
-- >>> jp = jackSymbolicPol 3 [2, 1] 'P'
-- >>> hallInnerProduct''' jp jp 5 == hallInnerProduct jp jp (constantRatioOfSprays 5)
hallInnerProduct''' :: 
  forall b. (Eq b, AlgField.C b, AlgMod.C (BaseRing b) b)
  => Spray b    -- ^ parametric spray
  -> Spray b    -- ^ parametric spray
  -> BaseRing b -- ^ parameter
  -> b 
hallInnerProduct''' = _hallInnerProduct psCombination (AlgMod.*>) 

-- | Hall inner product with parameter for parametric sprays. Same as 
-- @hallInnerProduct'''@ but with other constraints on the types. It is 
-- applicable to @SimpleParametricQSpray@ sprays, while @hallInnerProduct'''@ 
-- is not.
hallInnerProduct'''' :: 
  forall b. (Eq b, AlgRing.C b, AlgMod.C Rational b, AlgMod.C (BaseRing b) b)
  => Spray b    -- ^ parametric spray
  -> Spray b    -- ^ parametric spray
  -> BaseRing b -- ^ parameter
  -> b 
hallInnerProduct'''' = _hallInnerProduct psCombination' (AlgMod.*>) 

-- | the Hall inner product with symbolic parameter
_symbolicHallInnerProduct :: 
  (Eq a, AlgRing.C a) 
  => (Spray (Spray a) -> Spray (Spray a) -> Spray a -> Spray a) 
  -> Spray a -> Spray a -> Spray a
_symbolicHallInnerProduct func spray1 spray2 = func spray1' spray2' (lone 1)
  where
    spray1' = HM.map constantSpray spray1
    spray2' = HM.map constantSpray spray2

-- | Hall inner product with symbolic parameter. See README for some examples.
symbolicHallInnerProduct :: 
  (Eq a, AlgField.C a) => Spray a -> Spray a -> Spray a
symbolicHallInnerProduct =
  _symbolicHallInnerProduct 
    (
      _hallInnerProduct 
        (_psCombination (\spray_a r -> fromRational r *^ spray_a)) (^*^)
    ) 

-- | Hall inner product with symbolic parameter. Same as @symbolicHallInnerProduct@ 
-- but with other type constraints.
symbolicHallInnerProduct' :: 
  (Eq a, AlgMod.C Rational (Spray a), AlgRing.C a) 
  => Spray a -> Spray a -> Spray a
symbolicHallInnerProduct' =  _symbolicHallInnerProduct (hallInnerProduct')

-- | Hall inner product with symbolic parameter. Same as @symbolicHallInnerProduct@ 
-- but with other type constraints. It is applicable to @Spray Int@ sprays.
symbolicHallInnerProduct'' :: forall a. Real a => Spray a -> Spray a -> QSpray
symbolicHallInnerProduct'' spray1 spray2 = 
  _hallInnerProduct 
    (_psCombination (\qspray r -> r *^ qspray)) (^*^)
    qspray1' qspray2' (qlone 1)
  where
    asQSpray :: Spray a -> QSpray
    asQSpray = HM.map toRational
    qspray1' = HM.map constantSpray (asQSpray spray1)
    qspray2' = HM.map constantSpray (asQSpray spray2)

-- | Complete symmetric homogeneous polynomial
--
-- >>> putStrLn $ prettyQSpray (cshPolynomial 3 [2, 1])
-- 
cshPolynomial :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
cshPolynomial n lambda
  | n < 0                     = 
      error "cshPolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "cshPolynomial: invalid partition."
  | null lambda               = unitSpray
  | llambda > n               = zeroSpray
  | otherwise                 = productOfSprays (map cshPolynomialK lambda)
    where
      llambda = length lambda
      cshPolynomialK k = sumOfSprays msSprays
        where
          parts = partitions k
          msSprays = [msPolynomial n (fromPartition part) | part <- parts]

-- | power sum polynomial as a linear combination of 
-- complete symmetric homogeneous polynomials
pspInCSHbasis :: Partition -> Map Partition Rational
pspInCSHbasis mu = DM.fromList (zipWith f weights lambdas)
  where
    parts = partitions (sum mu)
    (weights, lambdas) = unzip $ filter ((/= 0) . fst) 
      [(eLambdaMu (fromPartition lambda) mu, fromPartition lambda) | lambda <- parts]
    f weight lambda = (lambda, weight)

-- | monomial symmetric polynomial as a linear combination of 
-- complete symmetric homogeneous polynomials
mspInCSHbasis :: Partition -> Map Partition Rational
mspInCSHbasis mu = sprayToMap (sumOfSprays sprays)
  where
    sprayToMap spray = 
      DM.fromList (HM.toList $ HM.mapKeys (DF.toList . exponents) spray) 
    comboToSpray combo = sumOfSprays 
      [ HM.singleton (Powers (S.fromList part) (length part)) c 
        | (part, c) <- DM.toList combo ]
    psAssocs = DM.toList (mspInPSbasis mu)
    sprays = 
      [c *^ comboToSpray (pspInCSHbasis lambda) | (lambda, c) <- psAssocs]

test :: Bool
test = msp == sumOfSprays cshs
  where
    mu = [3, 2, 1, 1]
    msp = msPolynomial 7 mu :: QSpray
    cshs = [c *^ cshPolynomial 7 lambda | (lambda, c) <- DM.toList (mspInCSHbasis mu)]


-- test'' :: (String, String)
-- test'' = (prettyParametricQSpray result, prettyParametricQSprayABCXYZ ["a"] ["b"] $ result)
--   where 
--     jsp = jackSymbolicPol' 3 [2, 1] 'P'
--     result = hallSymbolic'' jsp jsp
