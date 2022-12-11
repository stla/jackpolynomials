module Main where
import Math.Algebra.JackPol
import Math.Algebra.Jack
import Data.Ratio
import Math.Algebra.MultiPol
import Test.Tasty       (defaultMain, testGroup)
import Test.Tasty.HUnit (assertEqual, testCase)

main :: IO ()
main = defaultMain $
  testGroup "Tests"
  [ testCase "jackPol" $ do
      let jp = jackPol 2 [3, 1] (2%1)
      let v = evalPoly jp [1, 1]
      assertEqual ""
        v
        (48%1),

    testCase "jack" $ do
      assertEqual ""
        (jack [1, 1] [3, 1] (2%1))
        (48%1),

    testCase "schurPol" $ do
      let sp1 = schurPol 4 [4]
      let sp2 = schurPol 4 [3, 1]
      let sp3 = schurPol 4 [2, 2]
      let sp4 = schurPol 4 [2, 1, 1]
      let sp5 = schurPol 4 [1, 1, 1, 1]
      let v = evalPoly (sp1 ^+^ 3 *^ sp2 ^+^ 2*^ sp3 ^+^ 3*^ sp4 ^+^ sp5) [2, 2, 2, 2]
      assertEqual ""
        v
        4096,

    testCase "schur" $ do
      let sp1 = schur [1, 1, 1, 1] [4]
      let sp2 = schur [1, 1, 1, 1] [3, 1]
      let sp3 = schur [1, 1, 1, 1] [2, 2]
      let sp4 = schur [1, 1, 1, 1] [2, 1, 1]
      let sp5 = schur [1, 1, 1, 1] [1, 1, 1, 1]
      assertEqual ""
        (sp1 + 3*sp2 + 2*sp3 + 3*sp4 + sp5)
        256,

    testCase "zonalPol" $ do
      let zp1 = zonalPol 4 [3] :: Polynomial Rational
      let zp2 = zonalPol 4 [2, 1] :: Polynomial Rational
      let zp3 = zonalPol 4 [1, 1, 1] :: Polynomial Rational
      let v = evalPoly (zp1 ^+^ zp2 ^+^ zp3) [2, 2, 2, 2]
      assertEqual ""
        v
        512,

    testCase "zonal" $ do
      let zp1 = zonal [2%1, 2%1, 2%1, 2%1] [3]
      let zp2 = zonal [2%1, 2%1, 2%1, 2%1] [2, 1]
      let zp3 = zonal [2%1, 2%1, 2%1, 2%1] [1, 1, 1]
      assertEqual ""
        (zp1 + zp2 + zp3)
        512
  ]
