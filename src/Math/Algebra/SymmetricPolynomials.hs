{-|
Module      : Math.Algebra.SymmetricPolynomials
Description : More symmetric polynomials.
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

A Jack polynomial can have a very long expression in the canonical basis. 
A considerably shorter expression is obtained by writing the polynomial as 
a linear combination of the monomial symmetric polynomials instead, which is 
always possible since Jack polynomials are symmetric. This is the initial 
motivation of this module. But now it contains much more stuff dealing with 
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
  , esPolynomial
  -- * Decomposition of symmetric polynomials
  , msCombination
  , psCombination
  , psCombination'
  , psCombination''
  , cshCombination
  , cshCombination'
  , esCombination
  , esCombination'
  , schurCombination
  , schurCombination'
  , jackCombination
  , jackSymbolicCombination
  , jackSymbolicCombination'
  -- * Printing symmetric polynomials
  , prettySymmetricNumSpray
  , prettySymmetricQSpray
  , prettySymmetricQSpray'
  , prettySymmetricParametricQSpray
  , prettySymmetricSimpleParametricQSpray
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
  -- * Kostka-Foulkes polynomials
  , kostkaFoulkesPolynomial
  , kostkaFoulkesPolynomial'
  , skewKostkaFoulkesPolynomial
  , skewKostkaFoulkesPolynomial'
  -- * qt-Kostka polynomials
  , qtKostkaPolynomials
  , qtKostkaPolynomials'
  , qtSkewKostkaPolynomials
  , qtSkewKostkaPolynomials'
  -- * Hall-Littlewood polynomials
  , hallLittlewoodPolynomial
  , hallLittlewoodPolynomial'
  , transitionsSchurToHallLittlewood
  , skewHallLittlewoodPolynomial
  , skewHallLittlewoodPolynomial'
  -- * t-Schur polynomials
  , tSchurPolynomial
  , tSchurPolynomial'
  , tSkewSchurPolynomial
  , tSkewSchurPolynomial'
  -- * Macdonald polynomials
  , macdonaldPolynomial
  , macdonaldPolynomial'
  , skewMacdonaldPolynomial
  , skewMacdonaldPolynomial'
  , macdonaldJpolynomial
  , macdonaldJpolynomial'
  , skewMacdonaldJpolynomial
  , skewMacdonaldJpolynomial'
  , modifiedMacdonaldPolynomial
  , modifiedMacdonaldPolynomial'
  -- * Flagged Schur polynomials
  , flaggedSchurPol
  , flaggedSchurPol'
  , flaggedSkewSchurPol
  , flaggedSkewSchurPol'
  -- * Factorial Schur polynomials
  , factorialSchurPol
  , factorialSchurPol'
  , skewFactorialSchurPol
  , skewFactorialSchurPol'
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
import           Data.IntMap.Strict               ( 
                                                    IntMap
                                                  )
import qualified Data.IntMap.Strict               as IM
import           Data.List                        ( 
                                                    foldl1'
                                                  , foldl'
                                                  , nub
                                                  )
import           Data.List.Extra                  ( 
                                                    unsnoc
                                                  , allSame
                                                  , drop1
                                                  )
import           Data.Map.Merge.Strict            ( 
                                                    merge
                                                  , dropMissing
                                                  , zipWithMatched 
                                                  )
import           Data.Map.Strict                  ( 
                                                    Map
                                                  , unionsWith
                                                  , insert
                                                  )
import qualified Data.Map.Strict                  as DM
import           Data.Maybe                       ( fromJust )
import           Data.Ratio                       ( (%) )
import           Data.Sequence                    ( 
                                                    Seq (..)
                                                  , (|>)
                                                  )
import qualified Data.Sequence                    as S
import           Math.Algebra.Hspray              (
                                                    FunctionLike (..)
                                                  , (/^)
                                                  , (.^)
                                                  , Spray
                                                  , Powers (..)
                                                  , QSpray
                                                  , QSpray'
                                                  , ParametricSpray
                                                  , ParametricQSpray
                                                  , SimpleParametricSpray
                                                  , SimpleParametricQSpray
                                                  , isZeroSpray
                                                  , lone
                                                  , qlone
                                                  , lone'
                                                  , fromList
                                                  , getCoefficient
                                                  , getConstantTerm
                                                  , isConstant
                                                  , (%//%)
                                                  , (%/%)
                                                  , RatioOfSprays (..)
                                                  , RatioOfQSprays
                                                  , constantRatioOfSprays
                                                  , zeroRatioOfSprays
                                                  , prettyRatioOfQSpraysXYZ
                                                  , showNumSpray
                                                  , showQSpray
                                                  , showQSpray'
                                                  , showSpray
                                                  , prettyQSprayXYZ
                                                  , zeroSpray
                                                  , unitSpray
                                                  , productOfSprays
                                                  , sumOfSprays
                                                  , constantSpray
                                                  , allExponents
                                                  , asSimpleParametricSprayUnsafe
                                                  )
import           Math.Algebra.Jack.Internal       ( 
                                                    Partition
                                                  , _isPartition
                                                  , sprayToMap
                                                  , comboToSpray 
                                                  , _inverseKostkaMatrix
                                                  , _inverseSymbolicKostkaMatrix
                                                  , _kostkaFoulkesPolynomial
                                                  , _skewKostkaFoulkesPolynomial
                                                  , _hallLittlewoodPolynomialsInSchurBasis
                                                  , _transitionMatrixHallLittlewoodSchur
                                                  , skewHallLittlewoodP
                                                  , skewHallLittlewoodQ
                                                  , isSkewPartition
                                                  , flaggedSemiStandardYoungTableaux
                                                  , tableauWeight
                                                  , isIncreasing
                                                  , flaggedSkewTableaux
                                                  , skewTableauWeight
                                                  , macdonaldPolynomialP
                                                  , macdonaldPolynomialQ
                                                  , skewMacdonaldPolynomialP
                                                  , skewMacdonaldPolynomialQ
                                                  , chi_lambda_mu_rho
                                                  , clambda
                                                  , clambdamu
                                                  , macdonaldJinMSPbasis
                                                  , inverseKostkaNumbers
                                                  , skewSchurLRCoefficients
                                                  )
import           Math.Algebra.JackPol             ( 
                                                    schurPol
                                                  )
import           Math.Combinat.Compositions       ( compositions1 )
import           Math.Combinat.Partitions.Integer ( 
                                                    fromPartition
                                                  , toPartition
                                                  , toPartitionUnsafe
                                                  , partitions 
                                                  , partitionWidth
                                                  )
import           Math.Combinat.Partitions.Skew    ( 
                                                    mkSkewPartition
                                                  )
import           Math.Combinat.Permutations       ( permuteMultiset )
import           Math.Combinat.Tableaux           ( semiStandardYoungTableaux )
import           Math.Combinat.Tableaux.GelfandTsetlin ( kostkaNumbersWithGivenMu )
import           Math.Combinat.Tableaux.Skew      ( 
                                                    SkewTableau (..) 
                                                  , semiStandardSkewTableaux 
                                                  )

-- | monomial symmetric polynomial
msPolynomialUnsafe :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
msPolynomialUnsafe n lambda
  = fromList $ zip permutations coefficients
    where
      llambda      = length lambda
      permutations = permuteMultiset (lambda ++ replicate (n-llambda) 0)
      coefficients = repeat AlgRing.one

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
  | length lambda > n         = zeroSpray
  | otherwise                 = msPolynomialUnsafe n lambda

-- | Checks whether a spray defines a symmetric polynomial.
--
-- >>> -- note that the sum of two symmetric polynomials is not symmetric
-- >>> -- if they have different numbers of variables:
-- >>> spray = schurPol' 4 [2, 2] ^+^ schurPol' 3 [2, 1]
-- >>> isSymmetricSpray spray
isSymmetricSpray :: (AlgRing.C a, Eq a) => Spray a -> Bool
isSymmetricSpray spray = spray == spray' 
  where
    assocs = msCombination' spray
    n      = numberOfVariables spray
    spray' = foldl1' (^+^) 
      (
        map (\(lambda, x) -> x *^ msPolynomial n lambda) assocs
      )

-- | Symmetric polynomial as a linear combination of monomial symmetric polynomials.
msCombination :: AlgRing.C a => Spray a -> Map Partition a
msCombination spray = DM.fromList (msCombination' spray)

msCombination' :: AlgRing.C a => Spray a -> [(Partition, a)]
msCombination' spray = 
  map (\lambda -> let mu = DF.toList lambda in (mu, getCoefficient mu spray)) 
        lambdas
  where
    decreasing ys = and (S.zipWith (>=) ys (S.drop 1 ys))
    lambdas = filter decreasing (allExponents spray)

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
-- symmetric polynomials.
--
-- >>> putStrLn $ prettySymmetricParametricQSpray ["a"] $ jackSymbolicPol' 3 [3, 1, 1] 'J'
-- { [ 4*a^2 + 10*a + 6 ] }*M[3,1,1] + { [ 8*a + 12 ] }*M[2,2,1]
prettySymmetricParametricQSpray :: [String] -> ParametricQSpray -> String
prettySymmetricParametricQSpray letters spray = 
  showSpray (prettyRatioOfQSpraysXYZ letters) ("{ ", " }") 
            showSymmetricMonomials mspray
  where
    mspray = makeMSpray spray

-- | Prints a symmetric simple parametric spray as a linear combination of monomial 
-- symmetric polynomials.
prettySymmetricSimpleParametricQSpray :: 
  [String] -> SimpleParametricQSpray -> String
prettySymmetricSimpleParametricQSpray letters spray = 
  showSpray (prettyQSprayXYZ letters) ("(", ")") 
            showSymmetricMonomials mspray
  where
    mspray = makeMSpray spray

-- | Laplace-Beltrami operator on the space of homogeneous symmetric polynomials;
-- neither symmetry and homogeneity are checked.
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
--  | any (> n) lambda          = zeroSpray
--  | llambda > n               = zeroSpray
  | otherwise                 = productOfSprays sprays
    where
--      llambda = length lambda
      -- lists = [[(Powers (S.replicate (j-1) 0 |> k) j, AlgRing.one) | j <- [1 .. n]] | k <- lambda]
      -- expts = map (\j -> (j, S.replicate j 0)) [0 .. n-1]
      -- f k = map (\(j, sq) -> (Powers (sq |> k) j, AlgRing.one)) expts
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
          and (zipWith (>=) ys (drop1 ys))
--          and [ys !! i >= ys !! (i+1) | i <- [0 .. length ys - 2]]
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
mspInPSbasis kappa = 
  DM.filter (/= 0) (DM.fromDistinctAscList lambdas_and_weights) 
  where
    parts = partitions (sum kappa)
    lambdas_and_weights = 
      [let lambda = fromPartition part 
           weight = eLambdaMu kappa lambda in
         (lambda, weight / toRational (zlambda lambda))
        | part <- parts]

-- | the factor in the Hall inner product
zlambda :: Partition -> Int
zlambda lambda = p
  where
    parts = nub lambda
    table = [(j, sum [fromEnum (k == j) | k <- lambda]) | j <- parts]
    p =  
      product [factorial mj * part^mj | (part, mj) <- table]
    factorial n = product [2 .. n]

_symmPolyCombination :: 
    forall a b. (Eq a, AlgRing.C a) 
  => (Partition -> Map Partition b) 
  -> (a -> b -> a) 
  -> Spray a 
  -> Map Partition a
_symmPolyCombination mspInSymmPolyBasis func spray =
  if constantTerm == AlgAdd.zero 
    then symmPolyMap
    else insert [] constantTerm symmPolyMap
  where
    constantTerm = getConstantTerm spray
    msCombo = 
      msCombination (spray <+ (AlgAdd.negate constantTerm)) :: Map Partition a
    f :: Partition -> a -> Map Partition a
    f lambda coeff = 
      DM.map (func coeff) symmPolyCombo
      where
        symmPolyCombo = mspInSymmPolyBasis lambda :: Map Partition b
    symmPolyMap = 
      DM.filter (/= AlgAdd.zero) 
        (DM.foldlWithKey' 
          (\m lambda coeff -> DM.unionWith (AlgAdd.+) m (f lambda coeff)) 
            DM.empty msCombo)

-- | symmetric polynomial as a linear combination of power sum polynomials
_psCombination :: 
  (Eq a, AlgRing.C a) => (a -> Rational -> a) -> Spray a -> Map Partition a
_psCombination = _symmPolyCombination mspInPSbasis

-- | Symmetric polynomial as a linear combination of power sum polynomials. 
-- Symmetry is not checked.
psCombination :: 
  (Eq a, AlgField.C a) => Spray a -> Map Partition a
psCombination = _psCombination (\coef r -> coef AlgRing.* fromRational r)

-- | Symmetric polynomial as a linear combination of power sum polynomials. 
-- Same as @psCombination@ but with other constraints on the base ring of the spray.
psCombination' :: 
  (Eq a, AlgMod.C Rational a, AlgRing.C a) 
  => Spray a -> Map Partition a
psCombination' = _psCombination (flip (AlgMod.*>))

-- | Symmetric parametric spray as a linear combination of power sum polynomials. 
psCombination'' ::
    (FunctionLike b, Eq b, AlgRing.C b, AlgField.C (BaseRing b))
  => Spray b    -- ^ parametric spray
  -> Map Partition b  
psCombination'' = 
  _psCombination (\coef r -> fromRational r *^ coef)

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

-- | Hall inner product with Jack parameter, aka Jack-scalar product. It 
-- makes sense only for symmetric sprays, and the symmetry is not checked. 
hallInnerProduct :: 
  (Eq a, AlgField.C a)
  => Spray a   -- ^ spray
  -> Spray a   -- ^ spray
  -> a         -- ^ parameter
  -> a 
hallInnerProduct = _hallInnerProduct psCombination (AlgRing.*)

-- | Hall inner product with parameter. Same as @hallInnerProduct@ but 
-- with other constraints on the base ring of the sprays.
hallInnerProduct' :: 
  (Eq a, AlgMod.C Rational a, AlgRing.C a)
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
  (Eq b, AlgField.C b, AlgMod.C (BaseRing b) b)
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
  (Eq b, AlgRing.C b, AlgMod.C Rational b, AlgMod.C (BaseRing b) b)
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
-- x^3 + 2*x^2.y + 2*x^2.z + 2*x.y^2 + 3*x.y.z + 2*x.z^2 + y^3 + 2*y^2.z + 2*y.z^2 + z^3
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
          msSprays = 
            [msPolynomialUnsafe n (fromPartition part) 
              | part <- parts, partitionWidth part <= n]

-- | power sum polynomial as a linear combination of 
-- complete symmetric homogeneous polynomials
pspInCSHbasis :: Partition -> Map Partition Rational
pspInCSHbasis mu = 
  DM.filter (/= 0) (DM.fromDistinctAscList lambdas_and_weights)
  where
    parts = partitions (sum mu) 
    assoc kappa = 
      let kappa' = fromPartition kappa in (kappa', eLambdaMu kappa' mu)
    lambdas_and_weights = map assoc parts

-- | monomial symmetric polynomial as a linear combination of 
-- complete symmetric homogeneous polynomials
mspInCSHbasis :: Partition -> Map Partition Rational
mspInCSHbasis mu = sprayToMap (sumOfSprays sprays)
  where
    psAssocs = DM.toList (mspInPSbasis mu)
    sprays = 
      [c *^ comboToSpray (pspInCSHbasis lambda) | (lambda, c) <- psAssocs]

-- | symmetric polynomial as a linear combination of 
-- complete symmetric homogeneous polynomials
_cshCombination :: 
  (Eq a, AlgRing.C a) 
  => (a -> Rational -> a) -> Spray a -> Map Partition a
_cshCombination = _symmPolyCombination mspInCSHbasis

-- | Symmetric polynomial as a linear combination of complete symmetric 
-- homogeneous polynomials. Symmetry is not checked.
cshCombination :: 
  (Eq a, AlgField.C a) => Spray a -> Map Partition a
cshCombination = 
  _cshCombination (\coef r -> coef AlgRing.* fromRational r)

-- | Symmetric polynomial as a linear combination of complete symmetric homogeneous polynomials. 
-- Same as @cshCombination@ but with other constraints on the base ring of the spray.
cshCombination' :: 
  (Eq a, AlgMod.C Rational a, AlgRing.C a) 
  => Spray a -> Map Partition a
cshCombination' = _cshCombination (flip (AlgMod.*>))

-- | Elementary symmetric polynomial.
--
-- >>> putStrLn $ prettyQSpray (esPolynomial 3 [2, 1])
-- x^2.y + x^2.z + x.y^2 + 3*x.y.z + x.z^2 + y^2.z + y.z^2
esPolynomial :: (AlgRing.C a, Eq a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Spray a
esPolynomial n lambda
  | n < 0                     = 
      error "esPolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "esPolynomial: invalid partition."
  | null lambda               = unitSpray
  | l > n || any (>n) lambda  = zeroSpray
  | otherwise                 = productOfSprays (map esPolynomialK lambda)
    where
      l = length lambda
      esPolynomialK k = msPolynomialUnsafe n (replicate k 1)

-- | power sum polynomial as a linear combination of 
-- elementary symmetric polynomials
pspInESbasis :: Partition -> Map Partition Rational
pspInESbasis mu = 
  DM.filter (/= 0) (DM.fromDistinctAscList lambdas_and_weights)
  where
    wmu = sum mu
    parts = partitions wmu
    e = wmu - length mu
    e_is_even = even e
    negateIf = if e_is_even then id else negate 
    pair kappa = (kappa, negateIf (eLambdaMu kappa mu))
    lambdas_and_weights = 
      [pair (fromPartition part) | part <- parts]

-- | monomial symmetric polynomial as a linear combination of 
-- elementary symmetric polynomials
mspInESbasis :: Partition -> Map Partition Rational
mspInESbasis mu = sprayToMap (sumOfSprays sprays)
  where
    psAssocs = DM.toList (mspInPSbasis mu)
    sprays = 
      [c *^ comboToSpray (pspInESbasis lambda) | (lambda, c) <- psAssocs]

-- | symmetric polynomial as a linear combination of 
-- elementary symmetric polynomials
_esCombination :: 
  (Eq a, AlgRing.C a) 
  => (a -> Rational -> a) -> Spray a -> Map Partition a
_esCombination = _symmPolyCombination mspInESbasis

-- | Symmetric polynomial as a linear combination of elementary symmetric polynomials. 
-- Symmetry is not checked.
esCombination :: 
  (Eq a, AlgField.C a) => Spray a -> Map Partition a
esCombination = 
  _esCombination (\coef r -> coef AlgRing.* fromRational r)

-- | Symmetric polynomial as a linear combination of elementary symmetric polynomials. 
-- Same as @esCombination@ but with other constraints on the base ring of the spray.
esCombination' :: 
  (Eq a, AlgMod.C Rational a, AlgRing.C a) 
  => Spray a -> Map Partition a
esCombination' = _esCombination (flip (AlgMod.*>))

-- | complete symmetric homogeneous polynomial as a linear combination of 
-- Schur polynomials
cshInSchurBasis :: Int -> Partition -> Map Partition Int
cshInSchurBasis n mu = 
  DM.filterWithKey (\k _ -> length k <= n) 
                    (DM.mapKeys fromPartition kNumbers)
  where
    kNumbers = kostkaNumbersWithGivenMu (toPartitionUnsafe mu)

-- | symmetric polynomial as a linear combination of Schur polynomials
_schurCombination :: 
  forall a. (Eq a, AlgRing.C a) 
  => (a -> Rational -> a) -> Spray a -> Map Partition a
_schurCombination func spray =
  if constantTerm == AlgAdd.zero 
    then schurMap
    else insert [] constantTerm schurMap
  where
    constantTerm = getConstantTerm spray
    cshCombo = 
      _cshCombination func (spray <+ (AlgAdd.negate constantTerm))
    f :: Partition -> a -> Map Partition a
    f lambda coeff = 
      DM.map ((func coeff) . toRational) schurCombo
      where
        schurCombo = cshInSchurBasis (numberOfVariables spray) lambda 
    schurMap = 
      DM.filter (/= AlgAdd.zero) 
        (DM.foldlWithKey' 
          (\m lambda coeff -> DM.unionWith (AlgAdd.+) m (f lambda coeff)) 
            DM.empty cshCombo)

-- | Symmetric polynomial as a linear combination of Schur polynomials. 
-- Symmetry is not checked.
schurCombination :: 
  (Eq a, AlgField.C a) => Spray a -> Map Partition a
schurCombination = 
  _schurCombination (\coef r -> coef AlgRing.* fromRational r)

-- | Symmetric polynomial as a linear combination of Schur polynomials. 
-- Same as @schurCombination@ but with other constraints on the base ring of the spray.
schurCombination' :: 
  (Eq a, AlgMod.C Rational a, AlgRing.C a) 
  => Spray a -> Map Partition a
schurCombination' = _schurCombination (flip (AlgMod.*>))

-- | monomial symmetric polynomials in Jack polynomials basis
msPolynomialsInJackBasis :: 
  (Eq a, AlgField.C a)
  => a -> Char -> Int -> Int -> Map Partition (Map Partition a)
msPolynomialsInJackBasis alpha which n weight = 
   _inverseKostkaMatrix n weight alpha which

-- | monomial symmetric polynomials in Jack polynomials basis
msPolynomialsInJackSymbolicBasis :: 
  (Eq a, AlgField.C a) 
  => Char -> Int -> Int -> Map Partition (Map Partition (RatioOfSprays a))
msPolynomialsInJackSymbolicBasis which n weight = 
  _inverseSymbolicKostkaMatrix n weight which

-- | Symmetric polynomial as a linear combination of Jack polynomials with a 
-- given Jack parameter. Symmetry is not checked.
jackCombination :: 
  (Eq a, AlgField.C a)
  => a                      -- ^ Jack parameter
  -> Char                   -- ^ which Jack polynomials, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> Spray a                -- ^ spray representing a symmetric polynomial
  -> Map Partition a        -- ^ map representing the linear combination; a partition @lambda@ in the keys of this map corresponds to the term @coeff *^ jackPol' n lambda alpha which@, where @coeff@ is the value attached to this key and @n@ is the number of variables of the spray
jackCombination alpha which spray = 
  if not (which `elem` ['J', 'C', 'P', 'Q']) 
    then error "jackCombination: invalid character, must be 'J', 'C', 'P' or 'Q'."
    else
      _symmPolyCombination 
        (\lambda -> (combos IM.! (sum lambda)) DM.! lambda) 
          (AlgRing.*) spray
  where
    weights = filter (/= 0) (map DF.sum (allExponents spray))
    n = numberOfVariables spray
    combos = 
      IM.fromList 
        (zip weights (map (msPolynomialsInJackBasis alpha which n) weights))

-- | Symmetric polynomial as a linear combination of Jack polynomials with 
-- symbolic parameter. Symmetry is not checked.
jackSymbolicCombination :: 
     Char                   -- ^ which Jack polynomials, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> QSpray                 -- ^ spray representing a symmetric polynomial
  -> Map Partition RatioOfQSprays -- ^ map representing the linear combination; a partition @lambda@ in the keys of this map corresponds to the term @coeff *^ jackSymbolicPol' n lambda which@, where @coeff@ is the value attached to this key and @n@ is the number of variables of the spray
jackSymbolicCombination which qspray = 
  if not (which `elem` ['J', 'C', 'P', 'Q']) 
    then error "jackSymbolicCombination: invalid character, must be 'J', 'C', 'P' or 'Q'."
    else _symmPolyCombination 
      (\lambda -> (combos IM.! (sum lambda)) DM.! lambda) 
        (AlgRing.*) (HM.map constantRatioOfSprays qspray)
  where
    weights = filter (/= 0) (map DF.sum (allExponents qspray))
    n = numberOfVariables qspray
    combos = 
      IM.fromList 
      (zip weights (map (msPolynomialsInJackSymbolicBasis which n) weights))

-- | Symmetric parametric polynomial as a linear combination of Jack polynomials 
-- with symbolic parameter. 
-- Similar to @jackSymbolicCombination@ but for a parametric spray.
jackSymbolicCombination' :: 
  (Eq a, AlgField.C a)
  => Char                            -- ^ which Jack polynomials, @'J'@, @'C'@, @'P'@ or @'Q'@
  -> ParametricSpray a               -- ^ parametric spray representing a symmetric polynomial
  -> Map Partition (RatioOfSprays a) -- ^ map representing the linear combination; a partition @lambda@ in the keys of this map corresponds to the term @coeff *^ jackSymbolicPol' n lambda which@, where @coeff@ is the value attached to this key and @n@ is the number of variables of the spray
jackSymbolicCombination' which spray = 
  if not (which `elem` ['J', 'C', 'P', 'Q']) 
    then error "jackSymbolicCombination': invalid character, must be 'J', 'C', 'P' or 'Q'."
    else _symmPolyCombination 
      (\lambda -> (combos IM.! (sum lambda)) DM.! lambda) 
        (AlgRing.*) spray
  where
    weights = filter (/= 0) (map DF.sum (allExponents spray))
    n = numberOfVariables spray
    combos = 
      IM.fromList 
      (zip weights (map (msPolynomialsInJackSymbolicBasis which n) weights))

-- | Kostka-Foulkes polynomial of two given partitions. This is a univariate 
-- polynomial whose value at @1@ is the Kostka number of the two partitions.
kostkaFoulkesPolynomial :: 
  (Eq a, AlgRing.C a) => Partition -> Partition -> Spray a
kostkaFoulkesPolynomial lambda mu 
  | not (_isPartition lambda) = 
      error "kostkaFoulkesPolynomial: invalid partition."
  | not (_isPartition mu)     = 
      error "kostkaFoulkesPolynomial: invalid partition."
  | otherwise                 = 
      _kostkaFoulkesPolynomial lambda mu

-- | Kostka-Foulkes polynomial of two given partitions. This is a univariate 
-- polynomial whose value at @1@ is the Kostka number of the two partitions.
kostkaFoulkesPolynomial' :: Partition -> Partition -> QSpray
kostkaFoulkesPolynomial' = kostkaFoulkesPolynomial

-- | Skew Kostka-Foulkes polynomial. This is a univariate polynomial associated
-- to a skew partition and a partition, and its value at @1@ is the skew Kostka 
-- number associated to these partitions.
skewKostkaFoulkesPolynomial :: 
  (Eq a, AlgRing.C a) 
  => Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Partition -- ^ integer partition; the equality of the weight of this partition with the weight of the skew partition is a necessary condition to get a non-zero polynomial
  -> Spray a
skewKostkaFoulkesPolynomial lambda mu nu 
  | not (isSkewPartition lambda mu) =
     error "skewKostkaFoulkesPolynomial: invalid skew partition."
  | not (_isPartition nu) =
     error "skewKostkaFoulkesPolynomial: invalid partition."
  | otherwise = 
      _skewKostkaFoulkesPolynomial lambda mu nu

-- | Skew Kostka-Foulkes polynomial. This is a univariate polynomial associated
-- to a skew partition and a partition, and its value at @1@ is the skew Kostka 
-- number associated to these partitions.
skewKostkaFoulkesPolynomial' :: 
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Partition -- ^ integer partition; the equality of the weight of this partition with the weight of the skew partition is a necessary condition to get a non-zero polynomial
  -> QSpray
skewKostkaFoulkesPolynomial' = skewKostkaFoulkesPolynomial 

-- | qt-Kostka polynomials, aka Kostka-Macdonald polynomials. These are bivariate
-- symmetric polynomials usually denoted by \(K_{\lambda, \mu}(q,t)\) for two 
-- integer partitions \(\lambda\) and \(mu\), and \(q\) and \(t\) denote the 
-- variables. One obtains the Kostka-Foulkes polynomials by substituting \(q\) 
-- with \(0\). For a given partition \(\mu\), the function returns the polynomials
-- \(K_{\lambda, \mu}(q,t)\) for all partitions \(\lambda\) of the same weight as 
-- \(\mu\).
qtKostkaPolynomials :: 
  (Eq a, AlgField.C a) 
  => Partition 
  -> Map Partition (Spray a) 
qtKostkaPolynomials mu 
  | not (_isPartition mu) =
      error "qtKostkaPolynomials: invalid integer partition."
  | null mu =
      DM.singleton [] unitSpray
  | otherwise = 
      DM.map _numerator scs 
  where
    psCombo = macdonaldJinPSbasis mu
    t = lone' 2 
    den lambda = productOfSprays [unitSpray ^-^ t k | k <- lambda]
    msCombo lambda = 
      msCombination (psPolynomial (length lambda) lambda)
    ikn = inverseKostkaNumbers (sum mu)
    coeffs lambda = 
      let combo = msCombo lambda in
        DM.map 
          (\ikNumbers -> 
            DF.sum $ DM.intersectionWith (*) combo ikNumbers) 
          ikn
    scs = DM.foldlWithKey
      (\m lambda c -> 
        let den_lambda = den lambda in
          DM.unionWith (AlgAdd.+) m 
            (DM.map 
              (\ikNumber -> (ikNumber .^ c) %//% den_lambda) 
                (coeffs lambda)
            )
      )
      DM.empty psCombo

-- | qt-Kostka polynomials, aka Kostka-Macdonald polynomials. These are bivariate
-- symmetric polynomials usually denoted by \(K_{\lambda, \mu}(q,t)\) for two 
-- integer partitions \(\lambda\) and \(mu\), and \(q\) and \(t\) denote the 
-- variables. One obtains the Kostka-Foulkes polynomials by substituting \(q\) 
-- with \(0\). For a given partition \(\mu\), the function returns the polynomials
-- \(K_{\lambda, \mu}(q,t)\) for all partitions \(\lambda\) of the same weight as 
-- \(\mu\).
qtKostkaPolynomials' :: 
     Partition 
  -> Map Partition QSpray
qtKostkaPolynomials' = qtKostkaPolynomials

-- | Skew qt-Kostka polynomials. These are bivariate
-- symmetric polynomials usually denoted by \(K_{\lambda/\mu, \nu}(q,t)\) for two 
-- integer partitions \(\lambda\) and \(mu\) defining a skew partition, an 
-- integer partition \(\nu\), and \(q\) and \(t\) denote the 
-- variables. One obtains the skew Kostka-Foulkes polynomials by substituting \(q\) 
-- with \(0\). For given partitions \(\lambda\) and \(\mu\), the function returns 
-- the polynomials \(K_{\lambda/\mu, \nu}(q,t)\) for all partitions \(\nu\) of the 
-- same weight as the skew partition.
qtSkewKostkaPolynomials :: 
  (Eq a, AlgField.C a) 
  => Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Map Partition (Spray a)
qtSkewKostkaPolynomials lambda mu 
  | not (isSkewPartition lambda mu) =
      error "qtSkewKostkaPolynomials: invalid skew partition."
  | lambda == mu =
      DM.singleton [] unitSpray
  | otherwise = 
      DM.fromDistinctAscList (map spray nus)
  where
    lrCoeffs = skewSchurLRCoefficients lambda mu
    nus = partitions (sum lambda - sum mu)
    spray nu = 
      let nu' = fromPartition nu in
        (
          nu',
          foldl'
            (^+^) 
              zeroSpray
                (DM.intersectionWith (.^) lrCoeffs (qtKostkaPolynomials nu'))
        )

-- | Skew qt-Kostka polynomials. These are bivariate
-- symmetric polynomials usually denoted by \(K_{\lambda/\mu, \nu}(q,t)\) for two 
-- integer partitions \(\lambda\) and \(mu\) defining a skew partition, an 
-- integer partition \(\nu\), and \(q\) and \(t\) denote the 
-- variables. One obtains the skew Kostka-Foulkes polynomials by substituting \(q\) 
-- with \(0\). For given partitions \(\lambda\) and \(\mu\), the function returns 
-- the polynomials \(K_{\lambda/\mu, \nu}(q,t)\) for all partitions \(\nu\) of the 
-- same weight as the skew partition.
qtSkewKostkaPolynomials' :: 
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Map Partition QSpray
qtSkewKostkaPolynomials' = qtSkewKostkaPolynomials

-- | Hall-Littlewood polynomial of a given partition. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in a single parameter
-- usually denoted by \(t\). When substituting \(t\) with \(0\) in the 
-- Hall-Littlewood \(P\)-polynomials, one obtains the Schur polynomials.
hallLittlewoodPolynomial :: 
  (Eq a, AlgRing.C a) 
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Char      -- ^ which Hall-Littlewood polynomial, @'P'@ or @'Q'@
  -> SimpleParametricSpray a
hallLittlewoodPolynomial n lambda which 
  | n < 0 = 
      error "hallLittlewoodPolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "hallLittlewoodPolynomial: invalid partition."
  | not (which `elem` ['P', 'Q']) =
      error "hallLittlewoodPolynomial: last argument must be 'P' or 'Q'."
  | null lambda = 
      unitSpray
  | length lambda > n = 
      zeroSpray
  | otherwise = sumOfSprays sprays
    where
      coeffs = _hallLittlewoodPolynomialsInSchurBasis which lambda
      sprays = 
        DM.elems 
          (DM.mapWithKey 
            (\mu c -> c *^ (HM.map constantSpray (schurPol n mu))) coeffs)

-- | Hall-Littlewood polynomial of a given partition. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in a single parameter
-- usually denoted by \(t\). When substituting \(t\) with \(0\) in the 
-- Hall-Littlewood \(P\)-polynomials, one obtains the Schur polynomials.
hallLittlewoodPolynomial' :: 
     Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> Char      -- ^ which Hall-Littlewood polynomial, @'P'@ or @'Q'@
  -> SimpleParametricQSpray
hallLittlewoodPolynomial' = hallLittlewoodPolynomial

-- | Hall-Littlewood polynomials as linear combinations of Schur polynomials.
transitionsSchurToHallLittlewood :: 
     Int   -- ^ weight of the partitions of the Hall-Littlewood polynomials
  -> Char  -- ^ which Hall-Littlewood polynomials, @'P'@ or @'Q'@
  -> Map Partition (Map Partition (Spray Int))
transitionsSchurToHallLittlewood weight which 
  | weight <= 0                   = 
      error "transitionsHallLittlewoodToSchur: negative weight."
  | not (which `elem` ['P', 'Q']) =
      error "transitionsHallLittlewoodToSchur: the character must be 'P' or 'Q'."
  | otherwise                     = 
      _transitionMatrixHallLittlewoodSchur which weight

-- | Skew Hall-Littlewood polynomial of a given skew partition. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in a single parameter
-- usually denoted by \(t\). When substituting \(t\) with \(0\) in the skew
-- Hall-Littlewood \(P\)-polynomials, one obtains the skew Schur polynomials.
skewHallLittlewoodPolynomial :: (Eq a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Char      -- ^ which skew Hall-Littlewood polynomial, @'P'@ or @'Q'@
  -> SimpleParametricSpray a
skewHallLittlewoodPolynomial n lambda mu which 
  | n < 0 = 
      error "skewHallLittlewoodPolynomial: negative number of variables."
  | not (isSkewPartition lambda mu) = 
      error "skewHallLittlewoodPolynomial: invalid skew partition."
  | not (which `elem` ['P', 'Q']) =
      error "skewHallLittlewoodPolynomial: the character must be 'P' or 'Q'."
  | n == 0 = 
      if lambda == mu then unitSpray else zeroSpray
  | otherwise = 
      if which == 'P' 
        then skewHallLittlewoodP n (S.fromList lambda) (S.fromList mu)
        else skewHallLittlewoodQ n (S.fromList lambda) (S.fromList mu)
  
-- | Skew Hall-Littlewood polynomial of a given skew partition. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in a single parameter
-- usually denoted by \(t\). When substituting \(t\) with \(0\) in the skew
-- Hall-Littlewood \(P\)-polynomials, one obtains the skew Schur polynomials.
skewHallLittlewoodPolynomial' :: 
     Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Char      -- ^ which skew Hall-Littlewood polynomial, @'P'@ or @'Q'@
  -> SimpleParametricQSpray
skewHallLittlewoodPolynomial' = skewHallLittlewoodPolynomial

_tSkewSchurPolynomial ::
  (Eq a, AlgField.C a)
  => (Integer -> Integer -> a)
  -> Int
  -> Partition
  -> Partition
  -> SimpleParametricSpray a
_tSkewSchurPolynomial f n lambda mu = sumOfSprays sprays
  where
    w = sum lambda - sum mu
    rhos = partitions w
    t = lone' 1
    mapOfSprays = 
      IM.fromDistinctAscList (map (\r -> (r, unitSpray ^-^ t r)) [1 .. w])
    tPowerSumPol rho = 
      HM.map 
        (flip (*^) (productOfSprays (map ((IM.!) mapOfSprays) rho))) 
          (psPolynomial n rho)
    lambda' = S.fromList lambda
    mu' = S.fromList mu
    chi_lambda_mu_rhos = 
      [(rho', chi_lambda_mu_rho lambda' mu' (S.fromList rho')) 
        | rho <- rhos, let rho' = fromPartition rho]
    sprays = 
      [
        (f (toInteger c) (toInteger (zlambda rho)))
         AlgMod.*> tPowerSumPol rho
      | (rho, c) <- chi_lambda_mu_rhos, c /= 0
      ]

-- | t-Schur polynomial. This is a multivariate symmetric polynomial whose 
-- coefficients are polynomial in a single parameter usually denoted by \(t\).
-- One obtains the Schur polynomials by substituting \(t\) with \(0\). 
tSchurPolynomial ::
  (Eq a, AlgField.C a)
  => Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> SimpleParametricSpray a
tSchurPolynomial n lambda
  | n < 0 = 
      error "tSchurPolynomial: negative number of variables."
  | not (_isPartition lambda) =
      error "tSchurPolynomial: invalid partition."
  | otherwise =
      _tSkewSchurPolynomial 
        (\i j -> AlgRing.fromInteger i AlgField./ AlgRing.fromInteger j)
          n lambda []

-- | t-Schur polynomial. This is a multivariate symmetric polynomial whose 
-- coefficients are polynomial in a single parameter usually denoted by \(t\).
-- One obtains the Schur polynomials by substituting \(t\) with \(0\). 
tSchurPolynomial' ::
     Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> SimpleParametricQSpray
tSchurPolynomial' n lambda 
  | n < 0 = 
      error "tSchurPolynomial': negative number of variables."
  | not (_isPartition lambda) =
      error "tSchurPolynomial': invalid partition."
  | otherwise =
      _tSkewSchurPolynomial (%) n lambda []

-- | Skew t-Schur polynomial of a given skew partition. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in a single parameter
-- usually denoted by \(t\). One obtains the skew Schur polynomials by substituting 
-- \(t\) with \(0\). 
tSkewSchurPolynomial ::
  (Eq a, AlgField.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> SimpleParametricSpray a
tSkewSchurPolynomial n lambda mu
  | n < 0 = 
      error "tSkewSchurPolynomial: negative number of variables."
  | not (isSkewPartition lambda mu) = 
      error "tSkewSchurPolynomial: invalid skew partition."
  | otherwise =
      _tSkewSchurPolynomial 
        (\i j -> AlgRing.fromInteger i AlgField./ AlgRing.fromInteger j)
          n lambda mu

-- | Skew t-Schur polynomial of a given skew partition. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in a single parameter
-- usually denoted by \(t\). One obtains the skew Schur polynomials by substituting 
-- \(t\) with \(0\). 
tSkewSchurPolynomial' ::
     Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> SimpleParametricQSpray
tSkewSchurPolynomial' = _tSkewSchurPolynomial (%)

-- | Macdonald polynomial. This is a symmetric multivariate polynomial 
-- depending on two parameters usually denoted by \(q\) and \(t\).
-- Substituting \(q\) with \(0\) yields the Hall-Littlewood polynomials.
--
-- >>> macPoly = macdonaldPolynomial 3 [2, 1] 'P'
-- >>> putStrLn $ prettySymmetricParametricQSpray ["q", "t"] macPoly
-- { [ 1 ] }*M[2,1] + { [ 2*q.t^2 - q.t - q + t^2 + t - 2 ] %//% [ q.t^2 - 1 ] }*M[1,1,1]
macdonaldPolynomial :: (Eq a, AlgField.C a) 
  => Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> Char       -- ^ which Macdonald polynomial, @'P'@ or @'Q'@
  -> ParametricSpray a
macdonaldPolynomial n lambda which
  | n < 0 = 
      error "macdonaldPolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "macdonaldPolynomial: invalid partition."
  | not (which `elem` ['P', 'Q']) =
      error "macdonaldPolynomial: last argument must be 'P' or 'Q'."
  | null lambda = 
      unitSpray
  | length lambda > n = 
      zeroSpray
  | otherwise = 
      if which == 'P'
        then macdonaldPolynomialP n lambda
        else macdonaldPolynomialQ n lambda

-- | Macdonald polynomial. This is a symmetric multivariate polynomial 
-- depending on two parameters usually denoted by \(q\) and \(t\).
-- Substituting \(q\) with \(0\) yields the Hall-Littlewood polynomials.
macdonaldPolynomial' ::  
     Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> Char       -- ^ which Macdonald polynomial, @'P'@ or @'Q'@
  -> ParametricQSpray
macdonaldPolynomial' = macdonaldPolynomial

-- | Skew Macdonald polynomial of a given skew partition. This is a multivariate 
-- symmetric polynomial with two parameters usually denoted by \(q\) and \(t\).
-- Substituting \(q\) with \(0\) yields the skew Hall-Littlewood polynomials.
skewMacdonaldPolynomial :: (Eq a, AlgField.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Char      -- ^ which skew Macdonald polynomial, @'P'@ or @'Q'@
  -> ParametricSpray a
skewMacdonaldPolynomial n lambda mu which 
  | n < 0 = 
      error "skewMacdonaldPolynomial: negative number of variables."
  | not (isSkewPartition lambda mu) = 
      error "skewMacdonaldPolynomial: invalid skew partition."
  | not (which `elem` ['P', 'Q']) =
      error "skewMacdonaldPolynomial: the character must be 'P' or 'Q'."
  | n == 0 = 
      if lambda == mu then unitSpray else zeroSpray
  | otherwise = 
      if which == 'P' 
        then skewMacdonaldPolynomialP n lambda mu
        else skewMacdonaldPolynomialQ n lambda mu

-- | Skew Macdonald polynomial of a given skew partition. This is a multivariate 
-- symmetric polynomial with two parameters usually denoted by \(q\) and \(t\).
-- Substituting \(q\) with \(0\) yields the skew Hall-Littlewood polynomials.
skewMacdonaldPolynomial' :: 
     Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Char      -- ^ which skew Macdonald polynomial, @'P'@ or @'Q'@
  -> ParametricQSpray
skewMacdonaldPolynomial' = skewMacdonaldPolynomial

-- | Macdonald J-polynomial. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in two parameters.
macdonaldJpolynomial :: 
  forall a. (Eq a, AlgField.C a)
  => Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> SimpleParametricSpray a
macdonaldJpolynomial n lambda 
  | n < 0 = 
      error "macdonaldJpolynomial: negative number of variables."
  | not (_isPartition lambda) = 
      error "macdonaldJpolynomial: invalid partition."
  | null lambda = 
      unitSpray
  | length lambda > n = 
      zeroSpray
  | otherwise =
      asSimpleParametricSprayUnsafe $
        HM.map ((AlgMod.*>) (clambda (S.fromList lambda) :: Spray a)) 
          (macdonaldPolynomial n lambda 'P')

-- | Macdonald J-polynomial. This is a multivariate 
-- symmetric polynomial whose coefficients are polynomial in two parameters.
macdonaldJpolynomial' :: 
     Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> SimpleParametricQSpray
macdonaldJpolynomial' = macdonaldJpolynomial

-- | Skew Macdonald J-polynomial. This is a multivariate 
-- symmetric polynomial whose coefficients depend on two parameters.
skewMacdonaldJpolynomial :: 
  (Eq a, AlgField.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> ParametricSpray a
skewMacdonaldJpolynomial n lambda mu 
  | n < 0 = 
      error "skewMacdonaldJpolynomial: negative number of variables."
  | not (isSkewPartition lambda mu) = 
      error "skewMacdonaldJpolynomial: invalid skew partition."
  | n == 0 = 
      if lambda == mu then unitSpray else zeroSpray
  | otherwise = 
      clambdamu (S.fromList lambda) (S.fromList mu)  
        *^ skewMacdonaldPolynomial n lambda mu 'P'

-- | Skew Macdonald J-polynomial. This is a multivariate 
-- symmetric polynomial whose coefficients depend on two parameters.
skewMacdonaldJpolynomial' :: 
     Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> ParametricQSpray
skewMacdonaldJpolynomial' = skewMacdonaldJpolynomial

macdonaldJinPSbasis ::
  (Eq a, AlgField.C a) => Partition -> Map Partition (Spray a)
macdonaldJinPSbasis mu = 
  DM.filter (not . isZeroSpray) 
    (unionsWith (^+^) (DM.elems $ DM.mapWithKey combo_to_map macdonaldCombo))
  where
    macdonaldCombo = macdonaldJinMSPbasis mu
    combo_to_map lambda spray = 
      DM.map 
        (\r -> fromRational r *^ spray) 
          (mspInPSbasis lambda)

-- | Modified Macdonald polynomial. This is a multivariate symmetric polynomial
-- whose coefficients are polynomials in two parameters.
modifiedMacdonaldPolynomial :: 
     (Eq a, AlgField.C a) 
  => Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> SimpleParametricSpray a
modifiedMacdonaldPolynomial n mu 
  | n < 0 = 
      error "modifiedMacdonaldPolynomial: negative number of variables."
  | not (_isPartition mu) = 
      error "modifiedMacdonaldPolynomial: invalid partition."
  | null mu = 
      unitSpray
  | otherwise = 
      jp 
  where
    psCombo = macdonaldJinPSbasis mu
    q' = lone' 1
    t' = lone' 2
    num_and_den Empty = undefined
    num_and_den (e :<| Empty) = (q' e, unitSpray)
    num_and_den (e1 :<| (e2 :<| _)) = (q' e1, t' e2)
    rOS_from_term powers coeff = coeff *^ RatioOfSprays spray1 spray2
      where
        (spray1, spray2) = num_and_den (exponents powers)
    toROS spray = 
      HM.foldlWithKey' 
        (\ros powers coeff -> ros AlgAdd.+ rOS_from_term powers coeff) 
          zeroRatioOfSprays spray
    den lambda = productOfSprays [t' k ^-^ unitSpray | k <- lambda]
    nmu = sum (zipWith (*) [1 .. ] (drop 1 mu))
    jp = DM.foldlWithKey 
      (\spray lambda c -> 
          spray ^+^ 
            _numerator (toROS (t' (nmu + sum lambda) ^*^ c) %/% den lambda) 
              *^ psPolynomial n lambda)
      zeroSpray psCombo

-- | Modified Macdonald polynomial. This is a multivariate symmetric polynomial
-- whose coefficients are polynomials in two parameters.
modifiedMacdonaldPolynomial' :: 
     Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> SimpleParametricQSpray 
modifiedMacdonaldPolynomial' = modifiedMacdonaldPolynomial

-- | Flagged Schur polynomial. A flagged Schur polynomial is not symmetric 
-- in general.
flaggedSchurPol :: 
  (Eq a, AlgRing.C a) 
  => Partition -- ^ integer partition
  -> [Int]     -- ^ lower bounds
  -> [Int]     -- ^ upper bounds
  -> Spray a
flaggedSchurPol lambda as bs
  | not (_isPartition lambda) =
      error "flaggedSchurPol: invalid partition."
  | not (allSame [llambda, las, lbs]) = 
      error "flaggedSchurPol: the partition and the lists of lower bounds and upper bounds must have the same length."
  | llambda == 0 =
      unitSpray
  | not (isIncreasing as) = 
      error "flaggedSchurPol: the list of lower bounds is not increasing."
  | not (isIncreasing bs) = 
      error "flaggedSchurPol: the list of upper bounds is not increasing."
  | or (zipWith (>) as bs) = 
      error "flaggedSchurPol: lower bounds must be smaller than upper bounds."
  | otherwise = sumOfSprays sprays
    where
      llambda = length lambda
      las = length as
      lbs = length bs
      tableaux = flaggedSemiStandardYoungTableaux lambda as bs
      monomial tableau = 
        productOfSprays $ zipWith lone' [1 ..] (tableauWeight tableau)
      sprays = map monomial tableaux

-- | Flagged Schur polynomial. A flagged Schur polynomial is not symmetric 
-- in general.
flaggedSchurPol' :: 
     Partition -- ^ integer partition
  -> [Int]     -- ^ lower bounds
  -> [Int]     -- ^ upper bounds
  -> QSpray
flaggedSchurPol' = flaggedSchurPol

-- | Flagged skew Schur polynomial. A flagged skew Schur polynomial is not symmetric 
-- in general.
flaggedSkewSchurPol :: 
  (Eq a, AlgRing.C a) 
  => Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> [Int]     -- ^ lower bounds
  -> [Int]     -- ^ upper bounds
  -> Spray a
flaggedSkewSchurPol lambda mu as bs
  | not (isSkewPartition lambda mu) =
      error "flaggedSkewSchurPol: invalid skew partition."
  | not (allSame [llambda, las, lbs]) = 
      error "flaggedSkewSchurPol: the outer partition and the lists of lower bounds and upper bounds must have the same length."
  | not (isIncreasing as) = 
      error "flaggedSkewSchurPol: the list of lower bounds is not increasing."
  | not (isIncreasing bs) = 
      error "flaggedSkewSchurPol: the list of upper bounds is not increasing."
  | or (zipWith (>) as bs) = 
      error "flaggedSkewSchurPol: lower bounds must be smaller than upper bounds."
  | lambda == mu =
      unitSpray
  | otherwise = sumOfSprays sprays
    where
      llambda = length lambda
      las = length as
      lbs = length bs
      tableaux = flaggedSkewTableaux lambda mu as bs
      monomial tableau = 
        productOfSprays $ zipWith lone' [1 ..] (skewTableauWeight tableau)
      sprays = map monomial tableaux

-- | Flagged skew Schur polynomial. A flagged skew Schur polynomial is not symmetric 
-- in general.
flaggedSkewSchurPol' :: 
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> [Int]     -- ^ lower bounds
  -> [Int]     -- ^ upper bounds
  -> QSpray
flaggedSkewSchurPol' = flaggedSkewSchurPol

-- | Factorial Schur polynomial. See
-- [Kreiman's paper](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v15i1r84/pdf)
-- /Products of factorial Schur functions/ for the definition.
factorialSchurPol :: 
  (Eq a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ integer partition
  -> [a]       -- ^ the sequence denoted by \(y\) in the reference paper 
  -> Spray a
factorialSchurPol n lambda y 
  | n < 0 = 
      error "factorialSchurPol: negative number of variables." 
  | not (_isPartition lambda) =
      error "factorialSchurPol: invalid integer partition."
  | n == 0 = 
      if l == 0 then unitSpray else zeroSpray
  | otherwise = 
      sumOfSprays sprays
  where
    l = length lambda
    tableaux = semiStandardYoungTableaux n (toPartition lambda)
    lones = [lone i | i <- [1 .. n]]
    idx tableau i j = 
      let row = tableau !! (i-1) 
          a = row !! (j-1)
      in (a, a + j - i) 
    factor tableau i j = 
      let (a, k) = idx tableau i j in lones !! (a-1) <+ y !! (k-1)
    i_ = [1 .. l]
    ij_ = [(i, j) | i <- i_, j <- [1 .. lambda !! (i-1)]]
    factors tableau = [factor tableau i j | (i, j) <- ij_]
    spray tableau = productOfSprays (factors tableau)
    sprays = map spray tableaux

-- | Factorial Schur polynomial. See
-- [Kreiman's paper](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v15i1r84/pdf)
-- /Products of factorial Schur functions/ for the definition.
factorialSchurPol' :: 
     Int        -- ^ number of variables
  -> Partition  -- ^ integer partition
  -> [Rational] -- ^ the sequence denoted by \(y\) in the reference paper
  -> QSpray
factorialSchurPol' = factorialSchurPol

-- | Skew factorial Schur polynomial. See 
-- [Macdonald's paper](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/SLC/opapers/s28macdonald.pdf)
-- /Schur functions: theme and variations/, 6th variation, for the definition.
skewFactorialSchurPol :: 
  (Eq a, AlgRing.C a)
  => Int       -- ^ number of variables
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> IntMap a  -- ^ the sequence denoted by \(a\) in the reference paper
  -> Spray a
skewFactorialSchurPol n lambda mu y 
  | n < 0 = 
      error "skewFactorialSchurPol: negative number of variables." 
  | not (isSkewPartition lambda mu) =
      error "skewFactorialSchurPol: invalid skew integer partition."
  | n == 0 = 
      if lambda == mu then unitSpray else zeroSpray
  | otherwise = 
      sumOfSprays sprays
  where
    skewPartition = mkSkewPartition (toPartition lambda, toPartition mu)
    skewTableaux = semiStandardSkewTableaux n skewPartition
    getSkewTableau (SkewTableau x) = x
    lones = [lone i | i <- [1 .. n]]
    idx tableau i j = 
      let (offset, entries) = tableau !! (i-1) 
          a = entries !! (j-1)
      in (a, a + offset + j - i) 
    factor tableau i j = 
      let (a, k) = idx tableau i j in lones !! (a-1) <+ y IM.! k
    i_ = [1 .. length lambda]
    ij_ tableau = 
      [(i, j) | i <- i_, j <- [1 .. length (snd (tableau !! (i-1)))]]
    factors tableau = [factor tableau i j | (i, j) <- ij_ tableau]
    spray tableau = productOfSprays (factors (getSkewTableau tableau))
    sprays = map spray skewTableaux

-- | Skew factorial Schur polynomial. See 
-- [Macdonald's paper](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/SLC/opapers/s28macdonald.pdf)
-- /Schur functions: theme and variations/, 6th variation, for the definition.
skewFactorialSchurPol' :: 
     Int             -- ^ number of variables
  -> Partition       -- ^ outer partition of the skew partition
  -> Partition       -- ^ inner partition of the skew partition
  -> IntMap Rational -- ^ the sequence denoted by \(a\) in the reference paper
  -> QSpray
skewFactorialSchurPol' = skewFactorialSchurPol

