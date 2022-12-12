# jackpolynomials

Schur polynomials have applications in combinatorics and zonal polynomials have
applications in multivariate statistics. They are particular cases of
[Jack polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to evaluate these polynomials. It can also compute their symbolic form.

___

```haskell
import Math.Algebra.Jack
import Data.Ratio
jack [1, 1] [3, 1] (2%1)
-- 48 % 1
```

```haskell
import Math.Algebra.JackPol
import Data.Ratio
import Math.Algebra.Spray
jp = jackPol 2 [3, 1] (2%1)
prettySpray show "x" jp
-- "(18 % 1) * x^(1, 3) + (12 % 1) * x^(2, 2) + (18 % 1) * x^(3, 1)"
evalSpray jp [1, 1]
-- 48 % 1
```


## References

* I.G. Macdonald. *Symmetric Functions and Hall Polynomials*. Oxford Mathematical Monographs. The Clarendon Press Oxford University Press, New York, second edition, 1995.

* J. Demmel and P. Koev. *Accurate and efficient evaluation of Schur and Jack functions*. Mathematics of computations, vol. 75, n. 253, 223-229, 2005.

* Jack polynomials. <https://www.symmetricfunctions.com/jack.htm>.
