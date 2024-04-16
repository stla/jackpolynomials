module Main (main) where
import Math.Algebra.Hspray                      ( evalSymbolicSpray
                                                , Rational'
                                                )
import Math.Algebra.JackPol                     ( jackPol'
                                                )
import Math.Algebra.JackSymbolicPol             ( jackSymbolicPol' )
import Miniterion                               ( bench
                                                , bgroup
                                                , defaultMain
                                                , whnf )

n :: Int
n = 5

lambda :: [Int]
lambda = [4, 2, 2, 1]

alpha :: Rational
alpha = 2

alpha' :: Rational'
alpha' = 2


main :: IO ()
main = 
  defaultMain
    [ bgroup "Jack"
      [ bench "jackPol with the given alpha"       $ 
          whnf jackPol' n lambda alpha
      , bench "jackSymbolicPol"                    $ 
          whnf jackSymbolicPol' n lambda
      , bench "jackSymbolicPol evaluated at alpha" $ 
          whnf evalSymbolicSpray (jackSymbolicPol' n lambda) alpha'
      ]
    ]
