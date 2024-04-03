{-|
Module      : Math.Algebra.Jack.Auxiliary
Description : Some utilities for Jack polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

A Jack polynomial can have a very long expression which can be considerably 
reduced if the polynomial is written in the basis formed by the monomial 
symmetric polynomials instead. This is the motivation of this module.
-}

module Math.Algebra.Jack.Auxiliary
  ( msPolynomial
  , msCombination
  , prettySymmetricSpray
  ) where
import qualified Algebra.Ring                     as AlgRing
import           Data.Function                    ( on )
import           Data.List                        ( nub, sortBy )
import           Data.Map.Strict                  ( Map )
import qualified Data.Map.Strict                  as DM
import           Math.Algebra.Hspray              (
                                                    Spray
                                                  --, Powers
                                                  --, constantSpray
                                                  , fromList
                                                  , getCoefficient
                                                  --, getConstantTerm
                                                  , toList
                                                  , zeroSpray
                                                  )
import           Math.Algebra.Jack.Internal       ( Partition , _isPartition )
import           Math.Combinat.Permutations       ( permuteMultiset )
import           Math.Combinat.Partitions.Integer ( fromPartition, mkPartition )
import           Data.Text                        ( Text
                                                  , append
                                                  , cons
                                                  , intercalate
                                                  , pack
                                                  , snoc
                                                  , unpack
                                                  )

-- | Monomial symmetric polynomials
--
-- >>> putStrLn $ prettySpray' (msPolynomial 3 [2, 1])
-- (1) x1^2x2 + (1) x1^2x3 + (1) x1x2^2 + (1) x1x3^2 + (1) x2^2x3 + (1) x2x3^2
msPolynomial 
  :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
msPolynomial n lambda
  | n < 0                     = error "msPolynomial: negative number of variables."
  | not (_isPartition lambda) = error "msPolynomial: invalid partition"
  | llambda > n               = zeroSpray
  | otherwise                 = fromList $ zip permutations coefficients
    where
      llambda      = length lambda
      permutations = permuteMultiset (lambda ++ replicate (n-llambda) 0)
      coefficients = repeat AlgRing.one

-- | Symmetric polynomial as a linear combination of monomial symmetric polynomials; 
-- this function has been introduced mainly for usage in `prettySymmetricSpray`, but 
-- has been exported because it can be useful
msCombination :: AlgRing.C a => Spray a -> Map Partition a
msCombination spray = DM.fromList (msCombination' spray)

msCombination' :: AlgRing.C a => Spray a -> [(Partition, a)]
msCombination' spray = 
  map (\lambda -> (lambda, getCoefficient lambda spray)) lambdas
  where
    lambdas = nub $ map (fromPartition . mkPartition . fst) (toList spray)

-- | Prints a symmetric spray as a linear combination of monomial symmetric polynomials
--
-- >>> putStrLn $ prettySymmetricSpray $ jackPol' 3 [3,1,1] 2 'J'
-- (42 % 1) * M[3, 1, 1] + (28 % 1) * M[2, 2, 1]
prettySymmetricSpray :: (Show a, AlgRing.C a) => Spray a -> String
prettySymmetricSpray spray = unpack $ intercalate (pack " + ") termsText
  where
    assocs = msCombination' spray
    termsText     = 
      map termText (sortBy (flip compare `on` fst) assocs)
    termText assoc = append
      (snoc (snoc (cons '(' $ snoc coefText ')') ' ') '*')
      (prettyMonomial $ fst assoc)
      where
        coefText = pack $ show (snd assoc)
        prettyMonomial :: Partition -> Text -- [0,2,1] -> "M[0, 2, 1]"
        prettyMonomial lambda = append (pack " M") (cons '[' $ snoc text ']')
          where
            text = intercalate (pack ", ") (map (pack . show) lambda)
