module Math.Algebra.Jack.Auxiliary
  ( msPolynomial
  , msCombination
  ) where
import qualified Algebra.Ring                     as AlgRing
import           Data.List                        ( nub )
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

-- | Monomial symmetric polynomials
--
-- >>> putStrLn $ prettySpray' (msPolynomial 3 [2, 1])
-- (1) x1^2x2 + (1) x1^2x3 + ......
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

msCombination 
  :: AlgRing.C a
  => Spray a
  -> Map Partition a
msCombination spray = DM.fromList assocs
  where
    -- constantTerm = constantSpray (getConstantTerm spray) -- useless ?
    lambdas = nub $ map (fromPartition . mkPartition . fst) (toList spray)
    assocs  = map (\lambda -> (lambda, getCoefficient lambda spray)) lambdas

