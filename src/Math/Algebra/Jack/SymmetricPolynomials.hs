{-|
Module      : Math.Algebra.Jack.SymmetricPolynomials
Description : Some utilities for Jack polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

A Jack polynomial can have a very long expression which can be considerably 
reduced if the polynomial is written in the basis formed by the monomial 
symmetric polynomials instead. This is the motivation of this module.
-}

module Math.Algebra.Jack.SymmetricPolynomials
  ( isSymmetricSpray
  , msPolynomial
  , msCombination
  , prettySymmetricNumSpray
  , prettySymmetricQSpray
  , prettySymmetricQSpray'
  , prettySymmetricOneParameterQSpray
  ) where
import qualified Algebra.Ring                     as AlgRing
import qualified Data.Foldable                    as DF
import           Data.List                        ( foldl1', nub )
import           Data.Map.Strict                  ( Map )
import qualified Data.Map.Strict                  as DM
import           Data.Sequence                    ( Seq )
import           Math.Algebra.Hspray              (
                                                    (^+^)
                                                  , (*^)
                                                  , Spray
                                                  , QSpray
                                                  , QSpray'
                                                  , OneParameterQSpray
                                                  , fromList
                                                  , getCoefficient
                                                  , numberOfVariables
                                                  , prettyRatioOfQPolynomials
                                                  , showNumSpray
                                                  , showQSpray
                                                  , showQSpray'
                                                  , showSpray
                                                  , toList
                                                  , zeroSpray
                                                  )
import           Math.Algebra.Jack.Internal       ( Partition , _isPartition )
import           Math.Combinat.Permutations       ( permuteMultiset )
import           Math.Combinat.Partitions.Integer ( fromPartition, mkPartition )

-- | Monomial symmetric polynomials
--
-- >>> putStrLn $ prettySpray' (msPolynomial 3 [2, 1])
-- (1) x1^2.x2 + (1) x1^2.x3 + (1) x1.x2^2 + (1) x1.x3^2 + (1) x2^2.x3 + (1) x2.x3^2
msPolynomial :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
msPolynomial n lambda
  | n < 0                     = error "msPolynomial: negative number of variables."
  | not (_isPartition lambda) = error "msPolynomial: invalid partition."
  | llambda > n               = zeroSpray
  | otherwise                 = fromList $ zip permutations coefficients
    where
      llambda      = length lambda
      permutations = permuteMultiset (lambda ++ replicate (n-llambda) 0)
      coefficients = repeat AlgRing.one

-- | Checks whether a spray defines a symmetric polynomial; this is useless for 
-- Jack polynomials because they always are symmetric, but this module contains 
-- everything needed to build this function and it can be useful in another context
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
-- M[3, 1, 1] + M[2, 2, 1]
prettySymmetricNumSpray :: (Num a, Ord a, Show a, AlgRing.C a) => Spray a -> String
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

-- | Prints a symmetric symbolic spray as a linear combination of monomial symmetric polynomials
--
-- >>> putStrLn $ prettySymmetricOneParameterQSpray "a" $ jackOneParameterPol' 3 [3, 1, 1] 'J'
-- { 4*a^2 + 10*a + 6 }*M[3,1,1] + { 8*a + 12 }*M[2,2,1]
prettySymmetricOneParameterQSpray :: String -> OneParameterQSpray -> String
prettySymmetricOneParameterQSpray a spray = 
  showSpray (prettyRatioOfQPolynomials a) ("{ ", " }") 
            showSymmetricMonomials mspray
  where
    mspray = makeMSpray spray
