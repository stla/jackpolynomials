module Main where
import           Data.Ratio
import           Math.Algebra.Hspray
import           Math.Algebra.Jack
import           Math.Algebra.Jack.HypergeoPQ
import           Math.Algebra.JackPol
import           Math.HypergeoMatrix
import           Test.Tasty                     ( defaultMain
                                                , testGroup
                                                )
import           Test.Tasty.HUnit               ( assertEqual
                                                , testCase
                                                )

main :: IO ()
main = defaultMain $ testGroup

  "Tests"

  [ testCase "jackPol" $ do
    let jp = jackPol 2 [3, 1] (2 % 1)
        v  = evalSpray jp [1, 1]
    assertEqual "" v (48 % 1)

  , testCase "jack" $ do
    assertEqual "" (jack [1, 1] [3, 1] (2 % 1)) (48 % 1)

  , testCase "schurPol" $ do
    let sp1 = schurPol 4 [4]
        sp2 = schurPol 4 [3, 1]
        sp3 = schurPol 4 [2, 2]
        sp4 = schurPol 4 [2, 1, 1]
        sp5 = schurPol 4 [1, 1, 1, 1]
        v = evalSpray (sp1 ^+^ 3 *^ sp2 ^+^ 2 *^ sp3 ^+^ 3 *^ sp4 ^+^ sp5) [2, 2, 2, 2]
    assertEqual "" v 4096

  , testCase "schur" $ do
    let sp1 = schur [1, 1, 1, 1] [4]
        sp2 = schur [1, 1, 1, 1] [3, 1]
        sp3 = schur [1, 1, 1, 1] [2, 2]
        sp4 = schur [1, 1, 1, 1] [2, 1, 1]
        sp5 = schur [1, 1, 1, 1] [1, 1, 1, 1]
    assertEqual "" (sp1 + 3 * sp2 + 2 * sp3 + 3 * sp4 + sp5) 256

  , testCase "zonalPol" $ do
    let zp1 = zonalPol 4 [3]       :: Spray Rational
        zp2 = zonalPol 4 [2, 1]    :: Spray Rational
        zp3 = zonalPol 4 [1, 1, 1] :: Spray Rational
        v   = evalSpray (zp1 ^+^ zp2 ^+^ zp3) [2, 2, 2, 2]
    assertEqual "" v 512

  , testCase "zonal" $ do
    let zp1 = zonal [2 % 1, 2 % 1, 2 % 1, 2 % 1] [3]
        zp2 = zonal [2 % 1, 2 % 1, 2 % 1, 2 % 1] [2, 1]
        zp3 = zonal [2 % 1, 2 % 1, 2 % 1, 2 % 1] [1, 1, 1]
    assertEqual "" (zp1 + zp2 + zp3) 512

  , testCase "hypergeometric function" $ do
    let a  = [1 % 1, 2 % 1]
        b  = [3 % 1]
        x  = [1 % 5, 1 % 2]
        h1 = hypergeoPQ 10 a b x
    h2 <- hypergeomat 10 2 a b x
    assertEqual "" h1 h2
  ]
