# jackpolynomials

*Jack, zonal, Schur and skew Schur polynomials.*

<!-- badges: start -->
[![Stack-lts](https://github.com/stla/jackpolynomials/actions/workflows/Stack-lts.yml/badge.svg)](https://github.com/stla/jackpolynomials/actions/workflows/Stack-lts.yml)
[![Stack-nightly](https://github.com/stla/jackpolynomials/actions/workflows/Stack-nightly.yml/badge.svg)](https://github.com/stla/jackpolynomials/actions/workflows/Stack-nightly.yml)
<!-- badges: end -->

Schur polynomials have applications in combinatorics and zonal polynomials have
applications in multivariate statistics. They are particular cases of
[Jack polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to evaluate these polynomials. It can also compute their symbolic form.

___

```haskell
import Math.Algebra.Jack
jack' [1, 1] [3, 1] 2 'J'
-- 48 % 1
```

```haskell
import Math.Algebra.JackPol
import Math.Algebra.Hspray
jp = jackPol' 2 [3, 1] 2 'J'
putStrLn $ prettySpray' jp
-- (18 % 1) x1^3x2 + (12 % 1) x1^2x2^2 + (18 % 1) x1x2^3
evalSpray jp [1, 1]
-- 48 % 1
```

As of version `1.2.0.0`, it is possible to get Jack polynomials with a symbolic Jack parameter:

```haskell
import Math.Algebra.JackSymbolicPol
import Math.Algebra.Hspray
jp = jackSymbolicPol' 2 [3, 1] 'J'
putStrLn $ prettySymbolicQSpray "a" jp
-- ((2) + (4)a + (2)a^2)*x1^3x2 + ((4) + (4)a)*x1^2x2^2 + ((2) + (4)a + (2)a^2)*x1x2^3
putStrLn $ prettySpray' $ evalSymbolicSpray jp 2
-- (18 % 1) x1^3x2 + (12 % 1) x1^2x2^2 + (18 % 1) x1x2^3
```

From the definition of Jack polynomials, as well as from their implementation in this package, 
the coefficients of the Jack polynomials are fractions of polynomials in the Jack parameter. 
However, in the above example, one can see that the coefficients of the Jack polynomial `jp` 
are *polynomials* in the Jack parameter `a`. 
This fact actually is always true for the $J$-Jack polynomials (not for $P$ and $Q$). This is 
a consequence of the Knop & Sahi combinatorial formula.
But be aware that in spite of this fact, the coefficients of the polynomials returned by 
Haskell are *fractions* of polynomials (the type of these polynomials is `SymbolicSpray`, 
defined in the **hspray** package).


## References

* I.G. Macdonald. *Symmetric Functions and Hall Polynomials*. Oxford Mathematical Monographs. The Clarendon Press Oxford University Press, New York, second edition, 1995.

* J. Demmel and P. Koev. *Accurate and efficient evaluation of Schur and Jack functions*. Mathematics of computations, vol. 75, n. 253, 223-229, 2005.

* Jack polynomials. <https://www.symmetricfunctions.com/jack.htm>.
