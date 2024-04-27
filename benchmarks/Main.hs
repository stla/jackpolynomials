module Main (main) where
import Math.Algebra.Hspray                      ( Rational'
                                                , QSpray
                                                , QSpray'
                                                , ParametricQSpray
                                                , substituteParameters                                              
                                                )
import Math.Algebra.JackPol                     ( jackPol'
                                                )
import Math.Algebra.JackSymbolicPol             ( jackSymbolicPol' )
import Miniterion                               ( bench
                                                , bgroup
                                                , defaultMain
                                                , whnf )

nT :: Int
nT = 5

lambdaT :: [Int]
lambdaT = [4, 2, 2, 1]

alphaT :: Rational
alphaT = 2

alphaT' :: Rational'
alphaT' = 2

jP :: (Int, [Int], Rational) -> QSpray 
jP (n, lambda, alpha) = jackPol' n lambda alpha 'J'

jSP :: (Int, [Int]) -> ParametricQSpray
jSP (n, lambda) = jackSymbolicPol' n lambda 'J'

jSPeval :: (Int, [Int], Rational) -> QSpray
jSPeval (n, lambda, alpha) = substituteParameters (jackSymbolicPol' n lambda 'J') [alpha]

main :: IO ()
main = 
  defaultMain
    [ bgroup "Jack"
      [ bench "jackPol with the given alpha"       $ 
          whnf jP (nT, lambdaT, alphaT)
      , bench "jackSymbolicPol"                    $ 
          whnf jSP (nT, lambdaT)
      , bench "jackSymbolicPol evaluated at alpha" $ 
          whnf jSPeval (nT, lambdaT, alphaT)
      ]
    ]
