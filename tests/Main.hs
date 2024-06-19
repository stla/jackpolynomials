module Main ( main ) where
import qualified Algebra.Additive               as AlgAdd
import qualified Algebra.Module                 as AlgMod
import qualified Data.HashMap.Strict            as HM
import qualified Data.IntMap.Strict             as IM
import qualified Data.Map.Strict                as DM
import           Data.Matrix                    ( 
                                                  fromLists
                                                )
import Data.Ratio                               ( (%) )
import Math.Algebra.Hspray                      ( FunctionLike (..)
                                                , Spray, QSpray
                                                , SimpleParametricSpray
                                                , SimpleParametricQSpray
                                                , asSimpleParametricSpray
                                                , lone, qlone 
                                                , zeroSpray
                                                , unitSpray
                                                , evalSpray 
                                                , evalParametricSpray'
                                                , substituteParameters
                                                , canCoerceToSimpleParametricSpray
                                                , isHomogeneousSpray
                                                , asRatioOfSprays
                                                , constantRatioOfSprays
                                                , (%//%)
                                                , (/^)
                                                , sumOfSprays
                                                , productOfSprays
                                                , detLaplace
                                                , getConstantTerm
                                                )
import qualified Math.Algebra.Hspray            as Hspray
import Math.Algebra.Jack                        ( schur, skewSchur 
                                                , jack', zonal'
                                                , Partition )
import Math.Algebra.Jack.HypergeoPQ             ( hypergeoPQ )
import Math.Algebra.JackPol                     ( zonalPol, zonalPol', jackPol'
                                                , schurPol, schurPol', skewSchurPol' )
import Math.Algebra.JackSymbolicPol             ( jackSymbolicPol' )
import Math.Algebra.SymmetricPolynomials        ( isSymmetricSpray
                                                , prettySymmetricParametricQSpray
                                                , laplaceBeltrami
                                                , calogeroSutherland
                                                , hallInnerProduct
                                                , hallInnerProduct''
                                                , symbolicHallInnerProduct
                                                , symbolicHallInnerProduct''
                                                , msPolynomial
                                                , msCombination
                                                , psPolynomial 
                                                , psCombination
                                                , cshPolynomial
                                                , cshCombination
                                                , esPolynomial
                                                , esCombination
                                                , schurCombination
                                                , jackCombination
                                                , jackSymbolicCombination
                                                , jackSymbolicCombination'
                                                , kostkaNumbers
                                                , symbolicKostkaNumbers
                                                , kostkaFoulkesPolynomial
                                                , skewKostkaFoulkesPolynomial'
                                                , hallLittlewoodPolynomial
                                                , hallLittlewoodPolynomial'
                                                , skewHallLittlewoodPolynomial'
                                                , macdonaldPolynomial'
                                                , skewMacdonaldPolynomial'
                                                , flaggedSchurPol'
                                                , flaggedSkewSchurPol'
                                                , factorialSchurPol'
                                                , skewFactorialSchurPol'
                                                , tSchurPolynomial'
                                                , tSkewSchurPolynomial'
                                                , macdonaldJpolynomial'
                                                , skewMacdonaldJpolynomial'
                                                , qtKostkaPolynomials'
                                                , modifiedMacdonaldPolynomial'
                                                )
import Math.Combinat.Partitions.Integer         ( 
                                                  toPartition
                                                , fromPartition
                                                , mkPartition
                                                , partitions 
                                                , dualPartition
                                                )
import qualified Math.Combinat.Partitions.Integer as PI
import Math.Combinat.Tableaux.GelfandTsetlin    ( kostkaNumber )
import qualified Math.Combinat.Tableaux.GelfandTsetlin as GT
import Math.HypergeoMatrix                      ( hypergeomat )
import Test.Tasty                               ( defaultMain
                                                , testGroup
                                                )
import Test.Tasty.HUnit                         ( assertEqual
                                                , assertBool
                                                , testCase
                                                )

b_lambda_mu :: [Int] -> [Int] -> Int
b_lambda_mu lambda mu = sum $ DM.elems wholeMap -- zipWith (*) k1 k2 
  where
    parts = partitions (sum lambda)
    zeros = DM.fromList (zip parts (repeat 0))
    map1 = DM.union (GT.kostkaNumbersWithGivenMu (mkPartition lambda)) zeros
    map2 = DM.union (GT.kostkaNumbersWithGivenMu (mkPartition mu)) zeros
    wholeMap = DM.unionWithKey (\part kn1 _ -> kn1 * (map2 DM.! (dualPartition part))) map1 map2
    -- k1 = map ((flip kostkaNumber) (mkPartition lambda)) parts
    -- k2 = map (((flip kostkaNumber) (mkPartition mu)) . dualPartition) parts

a_lambda_mu :: [Int] -> [Int] -> Int
a_lambda_mu lambda mu = sum $ zipWith (*) k1 k2 
  where
    parts = partitions (sum lambda)
    k1 = map ((flip kostkaNumber) (mkPartition lambda)) parts
    k2 = map ((flip kostkaNumber) (mkPartition mu)) parts

main :: IO ()
main = defaultMain $ testGroup

  "Tests"

  [ 
  testCase "t-Schur polynomial" $ do
    let 
      tSchurPoly = tSchurPolynomial' 2 [2, 1]
      t = qlone 1
      x = lone 1 :: SimpleParametricQSpray
      y = lone 2 :: SimpleParametricQSpray
      expected = 
        (unitSpray ^-^ t) *^
          ((unitSpray ^-^ t)^**^2 *^ (x^**^2 ^*^ y ^+^ x ^*^ y^**^2)
            ^-^ t *^ (x^**^3 ^+^ y^**^3))
    assertEqual "" tSchurPoly expected

  , testCase "Skew t-Schur polynomial - branching rule" $ do
    let 
      lambda = [2, 2]
      tSchurPoly = tSchurPolynomial' 4 lambda
      ys = [lone 3, lone 4]
      expected = sumOfSprays
        [
          tSkewSchurPolynomial' 2 lambda mu 
            ^*^ changeVariables (tSchurPolynomial' 2 mu) ys 
          | mu <- [[], [1], [2], [1,1], [2,1], [2,2]]
        ]
    assertEqual "" tSchurPoly expected

  , testCase "Macdonald J-polynomial" $ do
    let
      n = 3
      macJpoly = macdonaldJpolynomial' n [2, 1]
      q = qlone 1
      t = qlone 2
      mus = [[3], [2, 1], [1, 1, 1]]
      qtKFpolys = [t, unitSpray ^+^ q ^*^ t, q]
      tSchurPolys = map ((HM.map (flip changeVariables [t])) . (tSchurPolynomial' n)) mus
      expected = sumOfSprays (zipWith (*^) qtKFpolys tSchurPolys)
    assertEqual "" macJpoly expected

  , testCase "Skew Macdonald J-polynomial at q=0" $ do
    let
      n = 3
      lambda = [3, 2]
      mu = [2, 1]
      skewHLpoly = skewHallLittlewoodPolynomial' n lambda mu 'Q'
      expected = asSimpleParametricSpray $
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing]))
          (skewMacdonaldJpolynomial' n lambda mu)
    assertEqual "" skewHLpoly expected

  , testCase "Macdonald polynomial branching rule" $ do
    let
      nx = 2
      ny = 2
      lambda = [2, 2]
      ys = [lone 3, lone 4]
      macJpoly = macdonaldJpolynomial' (nx + ny) lambda
      expected = asSimpleParametricSpray $
        sumOfSprays 
          [
            skewMacdonaldJpolynomial' nx lambda mu 
              ^*^ changeVariables (HM.map asRatioOfSprays $ macdonaldJpolynomial' ny mu) ys 
          | mu <- [[], [1], [2], [1,1], [2,1], [2,2]]
          ]
    assertEqual "" macJpoly expected

  , testCase "qt-Kostka polynomials" $ do
    let
      n = 4
      mu = [2, 1, 1]
      macJpoly = macdonaldJpolynomial' n mu
      qtKostkaPolys = qtKostkaPolynomials' mu
      expected = 
        sumOfSprays $ DM.elems $ DM.mapWithKey 
          (\lambda kp -> 
            kp *^ (HM.map (swapVariables (1, 2)) (tSchurPolynomial' n lambda)))
        qtKostkaPolys
    assertEqual "" macJpoly expected

  , testCase "Modified Macdonald polynomial" $ do
    let
      n = 4
      mu = [2, 1, 1]
      macHpoly = modifiedMacdonaldPolynomial' n mu
      macHpoly11 = substituteParameters macHpoly [1, 1] 
--      expected = productOfSprays [cshPolynomial n [k] | k <- mu] does not work
--      but with macHpoly = modifiedMacdonaldPolynomial' n [3, 1] (dual partition), this works
      expected = esPolynomial n [1] ^**^ (sum mu) 
    assertEqual "" macHpoly11 expected

  , testCase "Factorial Schur polynomial with y=[0 .. ] is Schur polynomial" $ do
    let 
      n = 4
      lambda = [3, 3, 2, 2]
      y = replicate (n + lambda !! 0 - 1) 0
      factorialSchurPoly = factorialSchurPol' n lambda y
      schurPoly = schurPol' n lambda
    assertEqual "" schurPoly factorialSchurPoly

  , testCase "Factorial Schur polynomial as determinant" $ do
    let 
      n = 3
      lambda = [3, 2, 2]
      y = [2, 6, 1, 2, 3]
      factorialSchurPoly = factorialSchurPol' n lambda y
      lones = [qlone i | i <- [1 .. n]]
      vandermonde = 
        productOfSprays [lones !! (i-1) ^-^ lones !! (j-1) 
                         | i <- [1 .. n-1], j <- [i+1 .. n]]
      x j k = productOfSprays [lones !! (j-1) <+ (y !! i) | i <- [0 .. k-1]]
      l = length lambda
      row i = [x i (lambda !! (j-1) + n - j) | j <- [1 .. l]]
      matrix = fromLists [row i | i <- [1 .. l]]
      det = detLaplace matrix
    assertEqual "" det (vandermonde ^*^ factorialSchurPoly)

  , testCase "Skew factorial Schur polynomial with y=0 is skew Schur polynomial" $ do
    let 
      n = 5
      lambda = [4, 3, 2, 2]
      mu = [2, 2]
      y = IM.fromList (zip [-2 .. 8] (repeat 0))
      skewFactorialSchurPoly = skewFactorialSchurPol' n lambda mu y
    assertEqual "" skewFactorialSchurPoly (skewSchurPol' n lambda mu)

  , testCase "Skew factorial Schur polynomial as determinant" $ do
    let 
      n = 3
      lambda = [3, 2, 2]
      mu = [2, 1]
      mu' = mu ++ [0]
      y = IM.fromList (zip [-2 ..] [2, 6, 1, 2, 3, 4, 5, 6])
      tau r = IM.mapKeys (subtract r) y
      skewFactorialSchurPoly = skewFactorialSchurPol' n lambda mu y
      kappa r = if r == 0 then [] else [r]
      h r a = if r < 0 then zeroSpray else factorialSchurPol' n (kappa r) a
      getSequence imap = [imap IM.! i | i <- [1 .. IM.size imap]]
      h' i j = h (lambda !! (i-1) - mu' !! (j-1) - i + j) 
                  (getSequence (tau (mu' !! (j-1) - j + 1)))
      l = length lambda
      row i = [h' i j | j <- [1 .. l]]
      matrix = fromLists [row i | i <- [1 .. l]]
      det = detLaplace matrix
    assertEqual "" det skewFactorialSchurPoly

  , testCase "Flagged Schur polynomial" $ do
    let 
      lambda = [5, 3, 2, 2]
      n = 5
      flaggedSchurPoly = flaggedSchurPol' lambda [1, 1, 1, 1] [n, n, n, n]
      schurPoly = schurPol' n lambda
    assertEqual "" flaggedSchurPoly schurPoly

  , testCase "Flagged skew Schur polynomial" $ do
    let 
      lambda = [5, 3, 2, 2]
      mu = [3, 1, 1]
      n = 5
      flaggedSkewSchurPoly = 
        flaggedSkewSchurPol' lambda mu [1, 1, 1, 1] [n, n, n, n]
      skewSchurPoly = skewSchurPol' n lambda mu
    assertEqual "" flaggedSkewSchurPoly skewSchurPoly

  , testCase "Jacobi-Trudi identity for flagged skew Schur polynomial" $ do
    let 
      lambda = [5, 3, 2, 2]
      mu = [3, 1, 1]
      as = [1, 1, 2, 4]
      bs = [2, 3, 4, 5]
      flaggedSkewSchurPoly = 
        flaggedSkewSchurPol' lambda mu as bs
      newVariables a b = map qlone [a .. b]
      h k a b
        | k < 0 = zeroSpray
        | k == 0 = changeVariables (cshPolynomial n []) variables
        | otherwise = changeVariables (cshPolynomial n [k]) variables
          where
            n = max 0 (b - a + 1)
            variables = newVariables a b
      l = length lambda
      mu' = mu ++ [0]
      row i = [h (lambda !! (i-1) - mu' !! (j-1) + j - i) (as !! (j-1)) (bs !! (i-1))| j <- [1 .. l]]
      matrix = fromLists [row i | i <- [1 .. l]]
      det = detLaplace matrix
    assertEqual "" det flaggedSkewSchurPoly

  , testCase "Jacobi-Trudi identity" $ do
    let 
      n = 5
      lambda = [3, 2, 1, 1]
      schurPoly = schurPol' n lambda
      h k 
        | k < 0 = zeroSpray
        | k == 0 = cshPolynomial n [] :: QSpray
        | otherwise = cshPolynomial n [k] :: QSpray
      row i = [h (lambda !! (i-1) + j - i) | j <- [1 .. 4]]
      matrix = fromLists [row i | i <- [1 .. 4]]
      det = detLaplace matrix
    assertEqual "" det schurPoly

  , testCase "jackSymbolicPol J" $ do
    let jp = jackSymbolicPol' 3 [3, 1] 'J'
        v  = evalParametricSpray' jp [2] [-3, 4, 5]
    assertEqual "" v 1488

  , testCase "jackSymbolicPol J has polynomial coefficients only" $ do
    let jp = jackSymbolicPol' 3 [3, 1] 'J'
    assertBool "" (canCoerceToSimpleParametricSpray jp)

  , testCase "jackSymbolicPol C" $ do
    let jp = jackSymbolicPol' 4 [3, 1] 'C'
        zp = zonalPol 4 [3, 1] :: Spray Rational
        p  = substituteParameters jp [2] 
    assertEqual "" zp p

  , testCase "jackSymbolicPol Q is symmetric" $ do
    let jp = jackSymbolicPol' 4 [3, 1] 'Q'
    assertBool "" (isSymmetricSpray jp)

  , testCase "jackSymbolicPol P is symmetric" $ do
    let jp = jackSymbolicPol' 5 [3, 2, 1] 'P'
    assertBool "" (isSymmetricSpray jp)

  , testCase "prettySymmetricParametricQSpray - jack J" $ do
    let jp = jackSymbolicPol' 3 [3, 1, 1] 'J'
    assertEqual "" 
      (prettySymmetricParametricQSpray ["a"] jp) 
      ("{ [ 4*a^2 + 10*a + 6 ] }*M[3,1,1] + { [ 8*a + 12 ] }*M[2,2,1]")

  , testCase "prettySymmetricParametricQSpray - jack C" $ do
    let jp = jackSymbolicPol' 3 [3, 1, 1] 'C'
    assertEqual "" 
      (prettySymmetricParametricQSpray ["a"] jp) 
      ("{ [ 20*a^2 ] %//% [ a^2 + (5/3)*a + (2/3) ] }*M[3,1,1] + { [ 40*a^2 ] %//% [ a^3 + (8/3)*a^2 + (7/3)*a + (2/3) ] }*M[2,2,1]")
 
  , testCase "jackPol" $ do
    let jp = jackPol' 2 [3, 1] (2 % 1) 'J'
        v  = evalSpray jp [1, 1]
    assertEqual "" v 48

  , testCase "jackPol is homogeneous" $ do
    let jp = jackPol' 4 [3, 1] (2 % 1) 'J'
    assertEqual "" (isHomogeneousSpray jp) (True, Just 4)

  , testCase "jackPol is symmetric (Groebner)" $ do
    let jp = jackPol' 3 [3, 2, 1] (2 % 1) 'J'
    assertBool "" (Hspray.isSymmetricSpray jp)

  , testCase "jack" $ do
    assertEqual "" (jack' [1, 1] [3, 1] (2 % 1) 'J') 48

  , testCase "Jack polynomial is eigenpolynomial for Laplace-Beltrami" $ do
    let
      alpha = 3 % 1
      lambda = [2, 2]
      b :: [Int] -> Rational
      b mu = toRational $ sum $ zipWith (*) mu [0 .. ]
      eigenvalue :: Int -> Rational -> [Int] -> Rational
      eigenvalue n a mu = 
        let mu' = fromPartition $ dualPartition (toPartition mu) in
          a * b mu' - b mu + toRational ((n-1) * sum mu)
      ev = eigenvalue 4 alpha lambda
      jp = jackPol' 4 lambda alpha 'J'
      jp' = laplaceBeltrami alpha jp
    assertEqual "" jp' (ev *^ jp)

  , testCase "Jack polynomial is eigenpolynomial for Calogero-Sutherland" $ do
    let
      eigenval :: Int -> Rational -> [Int] -> Rational
      eigenval n a mu = sum $ map 
        (\i -> let r = toRational (mu !! (i-1)) in 
                a/2 * r*r + ((toRational $ n + 1 - 2*i) / 2) * r) 
        [1 .. length mu]
      alpha = 3 % 4
      lambda = [3, 1]
      ev = eigenval 4 alpha lambda
      jp = jackPol' 4 lambda alpha 'J'
      jp' = calogeroSutherland alpha jp
    assertEqual "" jp' (ev *^ jp)

  , testCase "Jack P-polynomial for alpha=1 is Schur polynomial" $ do
    let 
      n = 5
      lambda = [5, 4, 3, 2, 1]
      jp = jackPol' n lambda 1 'P'
      sp = schurPol' n lambda
    assertEqual "" jp sp

  , testCase "schurPol" $ do
    let sp1 = schurPol 4 [4]
        sp2 = schurPol 4 [3, 1]
        sp3 = schurPol 4 [2, 2]
        sp4 = schurPol 4 [2, 1, 1]
        sp5 = schurPol 4 [1, 1, 1, 1] :: Spray Int
        v = evalSpray (sp1 ^+^ 3 *^ sp2 ^+^ 2 *^ sp3 ^+^ 3 *^ sp4 ^+^ sp5) 
                        [2, 2, 2, 2]
    assertEqual "" v 4096

  , testCase "schurPol is symmetric (Groebner)" $ do
    let sp = schurPol' 3 [3, 2, 1] 
    assertBool "" (Hspray.isSymmetricSpray sp)

  , testCase "schur" $ do
    let sp1 = schur [1, 1, 1, 1] [4]
        sp2 = schur [1, 1, 1, 1] [3, 1]
        sp3 = schur [1, 1, 1, 1] [2, 2]
        sp4 = schur [1, 1, 1, 1] [2, 1, 1]
        sp5 = schur [1, 1, 1, 1] [1, 1, 1, 1] :: Int
    assertEqual "" (sp1 + 3 * sp2 + 2 * sp3 + 3 * sp4 + sp5) 256

  , testCase "skewSchur" $ do
    let x = [2, 3, 4] :: [Int]
    assertEqual "" (skewSchur x [3, 2, 1] [1, 1]) 1890

  , testCase "skewSchurPol" $ do
    let x = lone 1 :: Spray Rational
        y = lone 2 :: Spray Rational
        z = lone 3 :: Spray Rational
        skp = skewSchurPol' 3 [2, 2, 1] [1, 1]
        p = x^**^2 ^*^ y  ^+^  x^**^2 ^*^ z  ^+^  x ^*^ y^**^2  
            ^+^  3 *^ (x ^*^ y ^*^ z)  ^+^  x ^*^ z^**^2  
            ^+^  y^**^2 ^*^ z  ^+^  y ^*^ z^**^2
    assertEqual "" skp p 

  , testCase "skewSchurPol is symmetric (Groebner)" $ do
    let skp = skewSchurPol' 3 [3, 2, 1] [1, 1]
    assertBool "" (Hspray.isSymmetricSpray skp)

  , testCase "zonalPol" $ do
    let zp1 = zonalPol' 4 [3]       
        zp2 = zonalPol' 4 [2, 1]    
        zp3 = zonalPol' 4 [1, 1, 1] 
        v   = evalSpray (zp1 ^+^ zp2 ^+^ zp3) [2, 2, 2, 2]
    assertEqual "" v 512

  , testCase "zonal" $ do
    let zp1 = zonal' [2 % 1, 2 % 1, 2 % 1, 2 % 1] [3]
        zp2 = zonal' [2 % 1, 2 % 1, 2 % 1, 2 % 1] [2, 1]
        zp3 = zonal' [2 % 1, 2 % 1, 2 % 1, 2 % 1] [1, 1, 1] 
    assertEqual "" (zp1 + zp2 + zp3) 512

  , testCase "hypergeometric function" $ do
    let a  = [1 % 1, 2 % 1]
        b  = [3 % 1]
        x  = [1 % 5, 1 % 2]
        h1 = hypergeoPQ 10 a b x :: Rational
    h2 <- hypergeomat 10 2 a b x
    assertEqual "" h1 h2

  , testCase "Hall inner product" $ do
    let
      alpha = 2
      poly1 = psPolynomial 4 [4]
      poly2 = psPolynomial 4 [3, 1]
      poly3 = psPolynomial 4 [2, 2]
      poly4 = psPolynomial 4 [2, 1, 1]
      poly5 = psPolynomial 4 [1, 1, 1, 1]
      h1 = hallInnerProduct poly1 poly1 alpha
      h2 = hallInnerProduct poly2 poly2 alpha
      h3 = hallInnerProduct poly3 poly3 alpha
      h4 = hallInnerProduct poly4 poly4 alpha
      h5 = hallInnerProduct poly5 poly5 alpha
      pow :: Rational -> Int -> Rational
      pow = (^)
    assertEqual ""
      (h1, h2, h3, h4, h5) 
      (
        4 * alpha
      , 3 * pow alpha 2
      , 8 * pow alpha 2
      , 4 * pow alpha 3
      , 24 * pow alpha 4
      )

  , testCase "Hall inner product of Jack P-polynomial and Jack Q-polynomial" $ do
    let
      jp1 = jackPol' 7 [4, 2, 1] 3 'P' 
      jp2 = jackPol' 7 [4, 2, 1] 3 'Q' 
      h = hallInnerProduct jp1 jp2 3
    assertEqual "" h 1

  , testCase "Hall inner product and b_lambda_mu" $ do
    let
      lambda = [4, 2, 1, 1]
      mu = [2, 2, 2, 2]
      h = cshPolynomial 8 lambda :: QSpray
      e = esPolynomial 8 mu :: QSpray
    assertEqual ""
      (hallInnerProduct h e 1) 
      (toRational $ b_lambda_mu lambda mu)

  , testCase "Hall inner product and a_lambda_mu" $ do
    let
      lambda = [4, 2, 2, 1]
      mu = [5, 4]
      hlambda = cshPolynomial 9 lambda :: QSpray
      hmu = cshPolynomial 9 mu :: QSpray
    assertEqual ""
      (hallInnerProduct hlambda hmu 1) 
      (toRational $ a_lambda_mu lambda mu)

  , testCase "Hall inner product of Schur polynomials" $ do
    let
      sp1 = schurPol 7 [4, 2, 1] :: Spray Int
      sp2 = schurPol 7 [2, 2, 2, 1] :: Spray Int
      h1 = hallInnerProduct'' sp1 sp1 1
      h2 = hallInnerProduct'' sp2 sp2 1
      h12 = hallInnerProduct'' sp1 sp2 1
    assertEqual "" (h1, h2, h12) (1, 1, 0)

  , testCase "Hall inner product with 'degenerate' symmetric polynomials" $ do
    let
      sp1 = schurPol' 3 [3,1]
      sp2 = schurPol' 3 [2,2]
      h1 = hallInnerProduct sp1 sp1 1
      h2 = hallInnerProduct sp2 sp2 1
      h12 = hallInnerProduct sp1 sp2 1
    assertEqual "" (h1, h2, h12) (10, 5, 6)

  , testCase "Symbolic Hall inner product with 'degenerate' symmetric polynomials" $ do
    let
      sp1 = schurPol' 3 [3,1]
      sp2 = schurPol' 3 [2,2]
      h1 = symbolicHallInnerProduct sp1 sp1
      h2 = symbolicHallInnerProduct sp2 sp2
      h12 = symbolicHallInnerProduct sp1 sp2
      alpha = qlone 1
    assertEqual "" 
      (h1, h2, h12) 
      (
        4*^alpha^**^3 ^+^ 5*^alpha^**^2 ^+^ alpha
      , alpha^**^3 ^+^ 3*^alpha^**^2 ^+^ alpha
      , 2*^alpha^**^3 ^+^ 3*^alpha^**^2 ^+^ alpha
      )

  , testCase "Symbolic Hall inner product" $ do
    let
      poly1 = psPolynomial 4 [4] :: QSpray
      poly2 = psPolynomial 4 [3, 1] :: QSpray
      poly3 = psPolynomial 4 [2, 2] :: QSpray
      poly4 = psPolynomial 4 [2, 1, 1] :: QSpray
      poly5 = psPolynomial 4 [1, 1, 1, 1] :: QSpray
      h1 = symbolicHallInnerProduct poly1 poly1
      h2 = symbolicHallInnerProduct poly2 poly2
      h3 = symbolicHallInnerProduct poly3 poly3
      h4 = symbolicHallInnerProduct poly4 poly4
      h5 = symbolicHallInnerProduct poly5 poly5
      alpha = qlone 1
    assertEqual ""
      (h1, h2, h3, h4, h5) 
      (
        4 *^ alpha
      , 3 *^ alpha^**^2
      , 8 *^ alpha^**^2
      , 4 *^ alpha^**^3
      , 24 *^ alpha^**^4
      )

  , testCase "Symbolic Hall inner product of Schur polynomials" $ do
    let
      sp1 = schurPol 7 [4, 2, 1] :: Spray Int
      sp2 = schurPol 7 [2, 2, 2, 1] :: Spray Int
      h1 = evaluateAt [1] (symbolicHallInnerProduct'' sp1 sp1)
      h2 = evaluateAt [1] (symbolicHallInnerProduct'' sp2 sp2)
      h12 = evaluateAt [1] (symbolicHallInnerProduct'' sp1 sp2)
    assertEqual "" (h1, h2, h12) (1, 1, 0)

  , testCase "Symbolic Hall inner product with parametric sprays" $ do
    let
      jp1 = jackSymbolicPol' 3 [1, 1, 1] 'P'
      jp2 = jackSymbolicPol' 3 [2, 1] 'P'
      jp3 = jackSymbolicPol' 3 [3] 'P'
      t = qlone 1
      t' = asRatioOfSprays t
      h1 = hallInnerProduct jp1 jp1 t'
      h2 = hallInnerProduct jp2 jp2 t'
      h3 = hallInnerProduct jp3 jp3 t'
    assertEqual ""
      (h1, h2, h3) 
      (
        asRatioOfSprays (t^**^3 /^ 6 ^+^ t^**^2 /^ 2 ^+^ t /^ 3)
      , (2*^t^**^3 ^+^ t^**^2) %//% (t <+ 2)
      , (3*^t^**^3) %//% ((t^**^2 ^+^ (3%2)*^t) <+ (1%2))
      )

  , testCase "Power sum polynomial and power sum combination" $ do
    let
      psPoly = 3*^psPolynomial 4 [2, 1, 1] ^-^ psPolynomial 4 [2, 1] :: QSpray
      psCombo = psCombination psPoly
    assertEqual ""
      psCombo 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "Complete symmetric homogeneous combination" $ do
    let
      cshPoly = 3*^cshPolynomial 4 [2, 1, 1] ^-^ cshPolynomial 4 [2, 1] :: QSpray
      cshCombo = cshCombination cshPoly
    assertEqual ""
      cshCombo 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "Elementary symmetric polynomials combination" $ do
    let
      esPoly = 3*^esPolynomial 4 [2, 1, 1] ^-^ esPolynomial 4 [2, 1] :: QSpray
      esCombo = esCombination esPoly
    assertEqual ""
      esCombo 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "Schur polynomials combination" $ do
    let
      poly = 3*^schurPol' 4 [2, 1, 1] ^-^ schurPol' 4 [2, 1]
      combo = schurCombination poly
    assertEqual ""
      combo 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "Schur polynomials combination of 'degenerate' symmetric polynomial" $ do
    let
      poly = psPolynomial 3 [4] :: QSpray
      schurCombo = schurCombination poly
      jackCombo = jackCombination 1 'P' poly
      expected = DM.fromList [([2, 1, 1], 1), ([3, 1], -1), ([4], 1)]
    assertEqual ""
      (schurCombo, jackCombo)
      (expected, expected)

  , testCase "Schur polynomials combination of a parametric spray" $ do
    let
      jpol = jackSymbolicPol' 4 [2, 2] 'J'
      schurCombo = schurCombination jpol
      alpha = qlone 1
      expected = [
          ([1, 1, 1, 1], asRatioOfSprays ((2*^alpha^**^2 ^-^ 6*^alpha) <+ 4)   )
        , ([2, 1, 1],    asRatioOfSprays (((-2)*^alpha^**^2 ^-^ 2*^alpha) <+ 4))
        , ([2, 2],       asRatioOfSprays ((2*^alpha^**^2 ^+^ 6*^alpha) <+ 4)   )
        ]
      jpol' = sumOfSprays $ map (\(lambda, c) -> c *^ schurPol 4 lambda) expected
    assertEqual ""
      (schurCombo, jpol) 
      (
        DM.fromList expected
      , jpol' 
      )

  , testCase "Jack J-polynomials combination" $ do
    let
      alpha = 3
      which = 'J'
      poly = 3*^jackPol' 4 [2, 1, 1] alpha which ^-^ jackPol' 4 [2, 1] alpha which
      combo = jackCombination alpha which poly
    assertEqual ""
      combo 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "Jack C-polynomials combination" $ do
    let
      alpha = 7
      which = 'C'
      poly = 3*^jackPol' 4 [2, 1, 1] alpha which ^-^ jackPol' 4 [2, 1] alpha which
      combo = jackCombination alpha which poly
    assertEqual ""
      combo 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "Jack Q-polynomials combination" $ do
    let
      which = 'Q'
      alpha = 4
      p = msPolynomial 5 [3, 1, 1] ^+^ psPolynomial 5 [3, 1] ^+^ 
          cshPolynomial 5 [2, 1] ^+^ esPolynomial 5 [2] ^+^ unitSpray :: QSpray
      sprays = [
          c *^ jackPol' 5 lambda alpha which 
          | (lambda, c) <- DM.toList (jackCombination alpha which p)
        ]
    assertEqual ""
      p (sumOfSprays sprays)

  , testCase "jackSymbolicCombination" $ do
    let
      alpha = 3
      which = 'J'
      poly = 3*^jackPol' 4 [2, 1, 1] alpha which ^-^ jackPol' 4 [2, 1] alpha which
      combo = jackSymbolicCombination which poly
      combo' = DM.filter (/= 0) (DM.map (evaluateAt [alpha]) combo)
    assertEqual ""
      combo' 
      (
        DM.fromList [([2, 1, 1], 3), ([2, 1], -1)]
      )

  , testCase "jackSymbolicCombination - 2" $ do
    let
      which = 'Q'
      p = 4 *^ msPolynomial 5 [3, 1, 1] ^+^ psPolynomial 5 [3, 1] ^-^ 
          5 *^ cshPolynomial 5 [2, 1] ^+^ unitSpray :: QSpray
      alpha = 7
      sprays = [
          (evaluateAt [alpha] c) *^ jackPol' 5 lambda alpha which 
          | (lambda, c) <- DM.toList (jackSymbolicCombination which p)
        ]
    assertEqual ""
      p (sumOfSprays sprays)

  , testCase "jackSymbolicCombination' (ParametricQSpray)" $ do
    let
      n = 4
      which = 'J'
      qspray = 7 *^ qlone 1
      poly = (3::Rational) AlgMod.*> jackSymbolicPol' n [2, 1, 1] which 
              ^+^ qspray AlgMod.*> jackSymbolicPol' n [2, 1] which
      combo = jackSymbolicCombination' which poly
      lambdas = DM.keys combo
      coeffs = DM.elems combo
    assertEqual ""
      (lambdas, coeffs)
      (
        [[2, 1], [2, 1, 1]]
      , [asRatioOfSprays qspray, constantRatioOfSprays 3]
      )

  , testCase "Kostka numbers" $ do
    let
      lambda = [4, 3, 1]
      kn1 = (kostkaNumbers (sum lambda) 1) DM.! lambda
      kn2 = DM.mapKeys fromPartition 
            (GT.kostkaNumbersWithGivenLambda (mkPartition lambda) :: DM.Map PI.Partition Rational)
    assertEqual "" kn1 kn2

  , testCase "Symbolic Kostka numbers" $ do
    let
      lambda = [4, 3, 1]
      kn1 = DM.map (evaluateAt [1]) (symbolicKostkaNumbers (sum lambda) DM.! lambda) :: DM.Map Partition Rational
      kn2 = DM.mapKeys fromPartition 
            (GT.kostkaNumbersWithGivenLambda (mkPartition lambda) :: DM.Map PI.Partition Rational)
    assertEqual "" kn1 kn2

  , testCase "Kostka-Foulkes polynomials" $ do
    let 
      lambda = [3, 1, 1]
      mu = [1, 1, 1, 1, 1]
      kfPoly = kostkaFoulkesPolynomial lambda mu :: Spray Int
      kNumber = kostkaNumber (toPartition lambda) (toPartition mu)
      kfPolyAt1 = evaluateAt [1] kfPoly
      t = lone 1 :: Spray Int
      expected = t^**^3 ^+^ t^**^4 ^+^ 2*^t^**^5 ^+^ t^**^6 ^+^ t^**^7
    assertEqual "" (kfPoly, kfPolyAt1) (expected, kNumber)

  , testCase "Skew Kostka-Foulkes polynomials" $ do
    let 
      lambda = [3, 3, 2, 1]
      mu = [1, 1, 1]
      n = sum lambda - sum mu
      nus = map fromPartition (partitions n)
      skewKFpolys = map (skewKostkaFoulkesPolynomial' lambda mu) nus
      hlPolys = map (\nu -> hallLittlewoodPolynomial' n nu 'P') nus
      combo = sumOfSprays $ zipWith (*^) skewKFpolys hlPolys
      comboIsConstant = all (== 0) (map numberOfVariables (HM.elems combo))
      combo' = HM.map getConstantTerm combo
      skewSchurPoly = skewSchurPol' n lambda mu
    assertEqual "" (comboIsConstant, combo') (True, skewSchurPoly)

  , testCase "Hall-Littlewood polynomial P" $ do
    let
      hlPoly = hallLittlewoodPolynomial 5 [2, 2, 1] 'P' :: SimpleParametricSpray Int
      msCombo = msCombination hlPoly
      t = lone 1 :: Spray Int
      expected = DM.fromList 
        [
          ([2, 2, 1], unitSpray)
        , ([2, 1, 1, 1], 2 +> (AlgAdd.negate (t ^+^ t^**^2)))
        , ([1, 1, 1, 1, 1], 5 +> ((-4)*^t ^-^ 4*^t^**^2 ^+^ t^**^3 ^+^ t^**^4 ^+^ t^**^5))
        ]
    assertEqual "" msCombo expected

  , testCase "Hall-Littlewood polynomial Q" $ do
    let
      hlQ2 = hallLittlewoodPolynomial 4 [2] 'Q' :: SimpleParametricSpray Int
      hlQ22 = hallLittlewoodPolynomial 4 [2, 2] 'Q'
      hlQ31 = hallLittlewoodPolynomial 4 [3, 1] 'Q'
      hlQ4 = hallLittlewoodPolynomial 4 [4] 'Q'
      spray = 1 +> (AlgAdd.negate (lone 1)) :: Spray Int
      expected = hlQ22 ^+^ spray *^ hlQ31 ^+^ spray *^ hlQ4
    assertEqual "" (hlQ2 ^**^ 2) expected

  , testCase "Skew Hall-Littlewood at t=0 is skew Schur polynomial" $ do
    let
      n = 3
      lambda = [3, 2, 1]
      mu = [1, 1]
      skewHLpoly = skewHallLittlewoodPolynomial' n lambda mu 'P'
      skewSchurPoly = skewSchurPol' n lambda mu
    assertEqual "" skewSchurPoly (substituteParameters skewHLpoly [0])

  , testCase "Skew Hall-Littlewood with mu=[] is Hall-Littlewood polynomial" $ do
    let
      n = 6
      lambda = [3, 2, 1]
      which = 'Q'
      skewHLpoly = skewHallLittlewoodPolynomial' n lambda [] which
      hlPoly = hallLittlewoodPolynomial' n lambda which
    assertEqual "" skewHLpoly hlPoly

  , testCase "Branching rule Hall-Littlewood P" $ do
    let
      lambda = [3, 1]
      mus = [[], [1], [2], [3], [1, 1], [2, 1], [3, 1]]
      nx = 2
      nz = 2
      which = 'P'
      hlLambda = hallLittlewoodPolynomial' (nx+nz) lambda which
      z = [lone 3, lone 4]
      terms = [skewHallLittlewoodPolynomial' nx lambda mu which ^*^ 
                changeVariables (hallLittlewoodPolynomial' nz mu which) z
                  | mu <- mus]
    assertEqual "" hlLambda (sumOfSprays terms)

  , testCase "Branching rule Hall-Littlewood Q" $ do
    let
      lambda = [3, 1]
      mus = [[], [1], [2], [3], [1, 1], [2, 1], [3, 1]]
      nx = 2
      nz = 2
      which = 'Q'
      hlLambda = hallLittlewoodPolynomial' (nx+nz) lambda which
      z = [lone 3, lone 4]
      terms = [skewHallLittlewoodPolynomial' nx lambda mu which ^*^ 
                changeVariables (hallLittlewoodPolynomial' nz mu which) z
                  | mu <- mus]
    assertEqual "" hlLambda (sumOfSprays terms)

  , testCase "Macdonald polynomial at q=0 is Hall-Littlewood" $ do
    let
      n = 3
      lambda = [2, 1]
      hlPolyP = hallLittlewoodPolynomial' n lambda 'P'
      hlPolyQ = hallLittlewoodPolynomial' n lambda 'Q'
      macPolyP = macdonaldPolynomial' n lambda 'P'
      macPolyQ = macdonaldPolynomial' n lambda 'Q'
      hlPolyP' = 
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing])) macPolyP
      hlPolyQ' = 
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing])) macPolyQ
    assertEqual "" 
      (hlPolyP, hlPolyQ)
      (asSimpleParametricSpray hlPolyP', asSimpleParametricSpray hlPolyQ')

  , testCase "Skew Macdonald polynomial at q=0 is skew Hall-Littlewood - 1" $ do
    let
      n = 3
      lambda = [3, 2]
      mu = [1, 1]
      shlPolyP = skewHallLittlewoodPolynomial' n lambda mu 'P'
      shlPolyQ = skewHallLittlewoodPolynomial' n lambda mu 'Q'
      smacPolyP = skewMacdonaldPolynomial' n lambda mu 'P'
      smacPolyQ = skewMacdonaldPolynomial' n lambda mu 'Q'
      shlPolyP' = 
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing])) smacPolyP
      shlPolyQ' = 
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing])) smacPolyQ
    assertEqual "" 
      (shlPolyP, shlPolyQ)
      (asSimpleParametricSpray shlPolyP', asSimpleParametricSpray shlPolyQ')

  , testCase "Skew Macdonald polynomial at q=0 is skew Hall-Littlewood - 2" $ do
    let
      n = 4
      lambda = [3, 2]
      mu = [1]
      shlPolyP = skewHallLittlewoodPolynomial' n lambda mu 'P'
      shlPolyQ = skewHallLittlewoodPolynomial' n lambda mu 'Q'
      smacPolyP = skewMacdonaldPolynomial' n lambda mu 'P'
      smacPolyQ = skewMacdonaldPolynomial' n lambda mu 'Q'
      shlPolyP' = 
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing])) smacPolyP
      shlPolyQ' = 
        HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing])) smacPolyQ
    assertEqual "" 
      (shlPolyP, shlPolyQ)
      (asSimpleParametricSpray shlPolyP', asSimpleParametricSpray shlPolyQ')

  ]
