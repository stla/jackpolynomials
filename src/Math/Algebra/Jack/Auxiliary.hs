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
-- (1) x1^2x2 + (1) x1^2x3 + ...... TODO
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

-- | Symmetric polynomial as a linear combination of monomial symmetric polynomials
msCombination 
  :: AlgRing.C a
  => Spray a
  -> Map Partition a
msCombination spray = DM.fromList (msCombination' spray)

msCombination' :: AlgRing.C a => Spray a -> [(Partition, a)]
msCombination' spray = 
  map (\lambda -> (lambda, getCoefficient lambda spray)) lambdas
  where
    lambdas = nub $ map (fromPartition . mkPartition . fst) (toList spray)
    -- constantTerm = constantSpray (getConstantTerm spray) -- useless ?

-- | prettyMonomial [0, 2, 1] = M[0, 2, 1]
prettyMonomial :: Partition -> Text
prettyMonomial lambda = append (pack " M") (cons '[' $ snoc string ']')
  where
    string = intercalate (pack ", ") (map (pack . show) lambda)

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
