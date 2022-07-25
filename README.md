# jackpolynomials

Schur polynomials have applications in combinatorics and zonal polynomials have
applications in multivariate statistics. They are particular cases of
[Jack polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to evaluate these polynomials. It can also compute their symbolic form.

```haskell
import Jack
import Data.Ratio
jackPol [1, 1] [3, 1] (2%1)
-- 48 % 1
```

```haskell
import JackPol
import Data.Ratio
import MultiPol
jp = jackPol 2 [3, 1] (2%1)
jp
-- (M (Monomial {coefficient = 18 % 1, powers = fromList [1,3]}) 
--  :+: 
--  M (Monomial {coefficient = 12 % 1, powers = fromList [2,2]})) 
--  :+: 
--  M (Monomial {coefficient = 18 % 1, powers = fromList [3,1]})
prettyPol show "x" jp
evalPoly jp [1, 1]
-- 48 % 1
```
