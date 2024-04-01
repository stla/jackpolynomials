module Main where
import Data.Ratio                               ( (%) )
import Math.Algebra.Hspray                      ( (^+^), (*^), (^*^), (^**^), Spray, lone
                                                , evalSpray, isSymmetricSpray )
import Math.Algebra.Jack                        ( jack, zonal, schur, skewSchur )
import Math.Algebra.Jack.HypergeoPQ             ( hypergeoPQ )
import Math.Algebra.JackPol                     ( zonalPol, jackPol, schurPol, skewSchurPol )
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

  [ testCase "jackPol" $ do
    let jp = jackPol 2 [3, 1] (2 % 1) "J" :: Spray Rational
        v  = evalSpray jp [1, 1]
    assertEqual "" v (48 % 1)

  , testCase "jackPol is symmetric" $ do
    let jp = jackPol 3 [3, 2, 1] (2 % 1) "J" :: Spray Rational
    assertBool "" (isSymmetricSpray jp)

  , testCase "jack" $ do
    assertEqual "" (jack [1, 1] [3, 1] (2 % 1)) (48 % 1 :: Rational)

  , testCase "schurPol" $ do
    let sp1 = schurPol 4 [4]
        sp2 = schurPol 4 [3, 1]
        sp3 = schurPol 4 [2, 2]
        sp4 = schurPol 4 [2, 1, 1]
        sp5 = schurPol 4 [1, 1, 1, 1] :: Spray Int
        v = evalSpray (sp1 ^+^ 3 *^ sp2 ^+^ 2 *^ sp3 ^+^ 3 *^ sp4 ^+^ sp5) [2, 2, 2, 2]
    assertEqual "" v 4096

  , testCase "schurPol is symmetric" $ do
    let sp = schurPol 3 [3, 2, 1] :: Spray Rational
    assertBool "" (isSymmetricSpray sp)

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
        skp = skewSchurPol 3 [2, 2, 1] [1, 1]
        p = x^**^2 ^*^ y  ^+^  x^**^2 ^*^ z  ^+^  x ^*^ y^**^2  ^+^  3 *^ (x ^*^ y ^*^ z) 
            ^+^  x ^*^ z^**^2  ^+^  y^**^2 ^*^ z  ^+^  y ^*^ z^**^2
    assertEqual "" skp p 

  , testCase "skewSchurPol is symmetric" $ do
    let skp = skewSchurPol 3 [3, 2, 1] [1, 1] :: Spray Rational
    assertBool "" (isSymmetricSpray skp)

  , testCase "zonalPol" $ do
    let zp1 = zonalPol 4 [3]       :: Spray Rational
        zp2 = zonalPol 4 [2, 1]    :: Spray Rational
        zp3 = zonalPol 4 [1, 1, 1] :: Spray Rational
        v   = evalSpray (zp1 ^+^ zp2 ^+^ zp3) [2, 2, 2, 2]
    assertEqual "" v 512

  , testCase "zonal" $ do
    let zp1 = zonal [2 % 1, 2 % 1, 2 % 1, 2 % 1] [3]
        zp2 = zonal [2 % 1, 2 % 1, 2 % 1, 2 % 1] [2, 1]
        zp3 = zonal [2 % 1, 2 % 1, 2 % 1, 2 % 1] [1, 1, 1] :: Rational
    assertEqual "" (zp1 + zp2 + zp3) 512

  , testCase "hypergeometric function" $ do
    let a  = [1 % 1, 2 % 1]
        b  = [3 % 1]
        x  = [1 % 5, 1 % 2]
        h1 = hypergeoPQ 10 a b x :: Rational
    h2 <- hypergeomat 10 2 a b x
    assertEqual "" h1 h2
  ]
