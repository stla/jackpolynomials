module Main ( main ) where
import Data.Ratio                               ( (%) )
import Math.Algebra.Hspray                      ( (^+^), (*^), (^*^), (^**^)
                                                , Spray, lone
                                                , evalSpray 
                                                , evalParametricSpray'
                                                , substituteParameters
                                                , canCoerceToSimpleParametricSpray
                                                , collinearSprays
                                                , isHomogeneousSpray
                                                )
import qualified Math.Algebra.Hspray            as Hspray
import Math.Algebra.Jack                        ( schur, skewSchur 
                                                , jack', zonal' )
import Math.Algebra.Jack.HypergeoPQ             ( hypergeoPQ )
import Math.Algebra.Jack.SymmetricPolynomials   ( isSymmetricSpray
                                                , prettySymmetricParametricQSpray
                                                , laplaceBeltrami
                                                , calogeroSutherland )
import Math.Algebra.JackPol                     ( zonalPol, zonalPol', jackPol'
                                                , schurPol, schurPol', skewSchurPol' )
import Math.Algebra.JackSymbolicPol             ( jackSymbolicPol' )
import Math.Combinat.Classes                    ( HasDuality (..) )
import Math.Combinat.Partitions.Integer         ( toPartition, fromPartition )
import Math.HypergeoMatrix                      ( hypergeomat )
import Test.Tasty                               ( defaultMain
                                                , testGroup
                                                )
import Test.Tasty.HUnit                         ( assertEqual
                                                , assertBool
                                                , testCase
                                                )

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
      eigenval n a mu = 
        sum $ map (\i -> let r = toRational (mu !! (i-1)) in a/2 * r*r + ((toRational $ n + 1 - 2*i) / 2) * r) [1 .. length mu]
      alpha = 3 % 4
      lambda = [3, 1]
      ev = eigenval 4 alpha lambda
      jp = jackPol' 4 lambda alpha 'J'
      jp' = calogeroSutherland alpha jp
    assertEqual "" jp' (ev *^ jp)

  , testCase "schurPol" $ do
    let sp1 = schurPol 4 [4]
        sp2 = schurPol 4 [3, 1]
        sp3 = schurPol 4 [2, 2]
        sp4 = schurPol 4 [2, 1, 1]
        sp5 = schurPol 4 [1, 1, 1, 1] :: Spray Int
        v = evalSpray (sp1 ^+^ 3 *^ sp2 ^+^ 2 *^ sp3 ^+^ 3 *^ sp4 ^+^ sp5) [2, 2, 2, 2]
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
        p = x^**^2 ^*^ y  ^+^  x^**^2 ^*^ z  ^+^  x ^*^ y^**^2  ^+^  3 *^ (x ^*^ y ^*^ z) 
            ^+^  x ^*^ z^**^2  ^+^  y^**^2 ^*^ z  ^+^  y ^*^ z^**^2
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
  ]
