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
motivation of this module. 
-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Math.Algebra.Jack.SymmetricPolynomials
  ( isSymmetricSpray
  , msPolynomial
  , msCombination
  , prettySymmetricNumSpray
  , prettySymmetricQSpray
  , prettySymmetricQSpray'
  , prettySymmetricParametricQSpray
  , laplaceBeltrami
  , calogeroSutherland
  , psPolynomial
  , psCombination
  , hallInnerProduct
  , psCombination'
  , hallInnerProduct'
  ) where
import           Prelude hiding ( fromIntegral, fromRational )
import qualified Algebra.Additive                 as AlgAdd
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
                                                  )
import           Math.Algebra.Jack.Internal       ( Partition , _isPartition )
import           Math.Combinat.Compositions       ( compositions1 )
import           Math.Combinat.Partitions.Integer ( 
                                                    fromPartition
                                                  , mkPartition
                                                  , partitions 
                                                  )
import           Math.Combinat.Permutations       ( permuteMultiset )

-- | Monomial symmetric polynomials
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

-- Laplace-Beltrami operator on the space of homogeneous symmetric polynomials;
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

-- Calogero-Sutherland operator on the space of homogeneous symmetric polynomials;
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

-- | Power sum polynomials
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
  | null lambda'              = unitSpray
  | llambda > n               = zeroSpray
  | otherwise                 = productOfSprays sprays
    where
      lambda' = fromPartition $ mkPartition lambda
      llambda      = length lambda'
      sprays = [HM.fromList $ [f i k | i <- [1 .. n]] | k <- lambda']
      f j k = (Powers expts j, AlgRing.one)
        where
          expts = S.replicate (j-1) 0 |> k

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
    ----
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
          [partitionSequences perm mu compo 
            | perm <- lambdaPerms, compo <- compos]
        xs = [eMuNus mu nus | nus <- sequencesOfPartitions]
    ----
    partitionSequences :: [Int] -> Partition -> [Int] -> [Partition]
    partitionSequences lambda mu compo = if test then nus else []
      where
        headOfCompo = fst $ fromJust (unsnoc compo)
        starts = scanl (+) 0 headOfCompo 
        ends   = zipWith (+) starts compo
        nus = [ 
                [ lambda !! k | k <- [starts !! i .. ends !! i - 1] ] 
                | i <- [0 .. length compo - 1]
              ]
        nuWeights = [sum nu | nu <- nus]
        decreasing xs = 
          and [xs !! i >= xs !! (i+1) | i <- [0 .. length xs - 2]]
        test = and (zipWith (==) mu nuWeights) && all decreasing nus
    ---- 
    eMuNus :: Partition -> [Partition] -> Rational
    eMuNus mu nus = product toMultiply
      where
        w :: Int -> Partition -> Rational
        w k nu = 
          let table = [sum [fromEnum (i == j) | i <- nu] | j <- nub nu] in
          (toInteger $ k * factorial (length nu - 1)) % 
            (toInteger $ product (map factorial table))
        factorial n = product [2 .. n]
        toMultiply = zipWith w mu nus

-- | the factor in the Hall inner product
zlambda :: Partition -> Int
zlambda lambda = p
  -- p AlgRing.* alpha AlgRing.^ (toInteger $ length lambda)
  where
    parts = nub lambda
    table = [sum [fromEnum (k == j) | k <- lambda] | j <- parts]
    p =  
      product [factorial mj * part^mj | (part, mj) <- zip parts table]
    factorial n = product [2 .. n]

-- | Symmetric polynomial as a linear combination of power sum polynomials
psCombination :: 
  forall a. (Eq a, AlgField.C a) => Spray a -> Map Partition a
psCombination spray =
  if constantTerm == AlgAdd.zero 
    then psMap
    else insert [] constantTerm psMap
  where
    constantTerm = getConstantTerm spray
    assocs = msCombination' (spray <+ (AlgAdd.negate constantTerm))
    f :: (Partition, a) -> [(Partition, a)] 
    f (lambda, coeff) = 
      map (second ((AlgRing.* coeff) . AlgField.fromRational)) (DM.toList psCombo)
      where
        psCombo = mspInPSbasis lambda :: Map Partition Rational
    psMap = DM.filter (/= AlgAdd.zero) 
            (unionsWith (AlgAdd.+) (map (DM.fromList . f) assocs))

-- | Symmetric polynomial as a linear combination of power sum polynomials
psCombination' :: 
  forall a. (Eq a, AlgMod.C Rational a, AlgRing.C a) => Spray a -> Map Partition a
psCombination' spray =
  if constantTerm == AlgAdd.zero 
    then psMap
    else insert [] constantTerm psMap
  where
    constantTerm = getConstantTerm spray
    assocs = msCombination' (spray <+ (AlgAdd.negate constantTerm))
    f :: (Partition, a) -> [(Partition, a)] 
    f (lambda, coeff) = 
      map (second (\r -> r AlgMod.*> coeff)) (DM.toList psCombo)
      where
        psCombo = mspInPSbasis lambda :: Map Partition Rational
    psMap = DM.filter (/= AlgAdd.zero) 
            (unionsWith (AlgAdd.+) (map (DM.fromList . f) assocs))

-- Hall inner product with parameter
hallInnerProduct :: 
  forall a. (Eq a, AlgField.C a)
  => Spray a   -- ^ spray
  -> Spray a   -- ^ spray
  -> a         -- ^ parameter
  -> a 
hallInnerProduct spray1 spray2 alpha = 
  AlgAdd.sum $ DM.elems
    (merge dropMissing dropMissing (zipWithMatched f) psCombo1 psCombo2)
  where
    psCombo1 = psCombination spray1 :: Map Partition a
    psCombo2 = psCombination spray2 :: Map Partition a
    zlambda' :: Partition -> a
    zlambda' lambda = fromIntegral (zlambda lambda) 
      AlgRing.* alpha AlgRing.^ (toInteger $ length lambda)
    f :: Partition -> a -> a -> a
    f lambda coeff1 coeff2 = 
      zlambda' lambda AlgRing.* (coeff1 AlgRing.* coeff2)

-- Hall inner product with parameter
hallInnerProduct' :: 
  forall a. (Eq a, AlgMod.C Rational a, AlgRing.C a)
  => Spray a   -- ^ spray
  -> Spray a   -- ^ spray
  -> a         -- ^ parameter
  -> a 
hallInnerProduct' spray1 spray2 alpha = 
  AlgAdd.sum $ DM.elems
    (merge dropMissing dropMissing (zipWithMatched f) psCombo1 psCombo2)
  where
    psCombo1 = psCombination' spray1 :: Map Partition a
    psCombo2 = psCombination' spray2 :: Map Partition a
    zlambda' :: Partition -> a
    zlambda' lambda = fromIntegral (zlambda lambda) 
      AlgRing.* alpha AlgRing.^ (toInteger $ length lambda)
    f :: Partition -> a -> a -> a
    f lambda coeff1 coeff2 = 
      zlambda' lambda AlgRing.* (coeff1 AlgRing.* coeff2)
