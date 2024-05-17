module Main ( main ) where
import qualified Data.Map.Strict                as DM
import Data.Ratio                               ( (%) )
import Math.Algebra.Hspray                      ( FunctionLike (..)
                                                , Spray, QSpray
                                                , lone, qlone 
                                                , evalSpray 
                                                , evalParametricSpray'
                                                , substituteParameters
                                                , canCoerceToSimpleParametricSpray
                                                , isHomogeneousSpray
                                                , asRatioOfSprays
                                                , (%//%)
                                                , (/^)
                                                )
import qualified Math.Algebra.Hspray            as Hspray
import Math.Algebra.Jack                        ( schur, skewSchur 
                                                , jack', zonal' )
import Math.Algebra.Jack.HypergeoPQ             ( hypergeoPQ )
import Math.Algebra.SymmetricPolynomials        ( isSymmetricSpray
                                                , prettySymmetricParametricQSpray
                                                , laplaceBeltrami
                                                , calogeroSutherland
                                                , hallInnerProduct
                                                , hallInnerProduct''
                                                , symbolicHallInnerProduct
                                                , symbolicHallInnerProduct''
                                                , psPolynomial 
                                                , psCombination
                                                , cshPolynomial
                                                , cshCombination
                                                , esPolynomial
                                                , esCombination
                                                , schurCombination
                                                )
import Math.Algebra.JackPol                     ( zonalPol, zonalPol', jackPol'
                                                , schurPol, schurPol', skewSchurPol' )
import Math.Algebra.JackSymbolicPol             ( jackSymbolicPol' )
import Math.Combinat.Classes                    ( HasDuality (..) )
import Math.Combinat.Partitions.Integer         ( 
                                                  toPartition
                                                , fromPartition
                                                , mkPartition
                                                , partitions 
                                                , dualPartition
                                                )
import Math.Combinat.Tableaux.GelfandTsetlin    ( kostkaNumber )
import Math.HypergeoMatrix                      ( hypergeomat )
import Test.Tasty                               ( defaultMain
                                                , testGroup
                                                )
import Test.Tasty.HUnit                         ( assertEqual
                                                , assertBool
                                                , testCase
                                                )

b_lambda_mu :: [Int] -> [Int] -> Int
b_lambda_mu lambda mu = sum $ zipWith (*) k1 k2 
  where
    parts = partitions (sum lambda)
    k1 = map ((flip kostkaNumber) (mkPartition lambda)) parts
    k2 = map (((flip kostkaNumber) (mkPartition mu)) . dualPartition) parts

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
  testCase "jackSymbolicPol J" $ do
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
        let mu' = fromPartition $ dual (toPartition mu) in
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

  ]
