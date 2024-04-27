# jackpolynomials

***Jack, zonal, Schur and skew Schur polynomials.***

<!-- badges: start -->
[![Stack-lts](https://github.com/stla/jackpolynomials/actions/workflows/Stack-lts.yml/badge.svg)](https://github.com/stla/jackpolynomials/actions/workflows/Stack-lts.yml)
[![Stack-nightly](https://github.com/stla/jackpolynomials/actions/workflows/Stack-nightly.yml/badge.svg)](https://github.com/stla/jackpolynomials/actions/workflows/Stack-nightly.yml)
<!-- badges: end -->

Schur polynomials have applications in combinatorics and zonal polynomials have
applications in multivariate statistics. They are particular cases of
[Jack polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to evaluate these polynomials and to compute them in symbolic form.

___

Evaluation of the Jack polynomial with parameter `2` associated to the integer 
partition `[3, 1]`, at `x1 = 1` and `x2 = 1`:

```haskell
import Math.Algebra.Jack
jack' [1, 1] [3, 1] 2 'J'
-- 48 % 1
```

The non-evaluated Jack polynomial:

```haskell
import Math.Algebra.JackPol
import Math.Algebra.Hspray
jp = jackPol' 2 [3, 1] 2 'J'
putStrLn $ prettyQSpray jp
-- 18*x^3.y + 12*x^2.y^2 + 18*x.y^3
evalSpray jp [1, 1]
-- 48 % 1
```

The first argument, here `2`, is the number of variables of the polynomial.


### Symbolic Jack parameter

As of version `1.2.0.0`, it is possible to get Jack polynomials with a 
symbolic Jack parameter:

```haskell
import Math.Algebra.JackSymbolicPol
import Math.Algebra.Hspray
jp = jackSymbolicPol' 2 [3, 1] 'J'
putStrLn $ prettyParametricQSpray jp
-- { [ 2*a^2 + 4*a + 2 ] }*X^3.Y + { [ 4*a + 4 ] }*X^2.Y^2 + { [ 2*a^2 + 4*a + 2 ] }*X.Y^3
putStrLn $ prettyQSpray' $ substituteParameters jp [2]
-- 18*x^3.y + 12*x^2.y^2 + 18*x.y^3
```

This is possible thanks to the **hspray** package which provides the type 
`ParametricSpray`. An object of this type represents a multivariate polynomial 
whose coefficients depend on some parameters which are symbolically treated. 
The type of the Jack polynomial returned by the `jackSymbolicPol` function is 
`ParametricSpray a`, and it is `ParametricQSpray` for the `jackSymbolicPol'` 
function. The type `ParametricQSpray` is an alias of `ParametricSpray Rational`.

From the definition of Jack polynomials, as well as from their implementation 
in this package, the coefficients of the Jack polynomials are 
*fractions of polynomials* in the Jack parameter. However, in the above 
example, one can see that the coefficients of the Jack polynomial `jp` are 
*polynomials* in the Jack parameter `a`. This fact actually is always true for 
the $J$-Jack polynomials (not for $C$, $P$ and $Q$). This is a consequence of 
the Knop & Sahi combinatorial formula. But be aware that in spite of this fact, 
the coefficients of the polynomials returned by Haskell are *fractions* of 
polynomials, in the sense that this is the nature of the `ParametricSpray` 
objects. 

Note that if you use the function `jackSymbolicPol` to get a 
`ParametricSpray Double` object in the output, it is not guaranted that you 
will visually get some polynomials in the Jack parameter for the coefficients, 
because the arithmetic operations are not exact with the `Double` type


### Showing symmetric polynomials

As of version 1.2.1.0, there is a module providing some functions to print a 
symmetric polynomial as a linear combination of the monomial symmetric 
polynomials. This can considerably shorten the expression of a symmetric 
polynomial as compared to its expression in the canonical basis, and the 
motivation to add this module to the package is that any Jack polynomial is 
a symmetric polynomial. Here is an example:

```haskell
import Math.Algebra.JackPol
import Math.Algebra.Jack.SymmetricPolynomials
jp = jackPol' 3 [3, 1, 1] 2 'J'
putStrLn $ prettySymmetricQSpray jp
-- 42*M[3,1,1] + 28*M[2,2,1]
```

And another example, with a symbolic Jack polynomial:

```haskell
import Math.Algebra.JackSymbolicPol
import Math.Algebra.Jack.SymmetricPolynomials
jp = jackSymbolicPol' 3 [3, 1, 1] 'J'
putStrLn $ prettySymmetricParametricQSpray ["a"] jp
-- { [ 4*a^2 + 10*a + 6 ] }*M[3,1,1] + { [ 8*a + 12 ] }*M[2,2,1]
```

Of course you can use these functions for other polynomials, but carefully: 
they do not check the symmetry. This new module provides the function 
`isSymmetricSpray` to check the symmetry of a polynomial, much more efficient 
than the function with the same name in the **hspray** package.


## References

* I.G. Macdonald. *Symmetric Functions and Hall Polynomials*. Oxford Mathematical Monographs. The Clarendon Press Oxford University Press, New York, second edition, 1995.

* J. Demmel and P. Koev. *Accurate and efficient evaluation of Schur and Jack functions*. Mathematics of computations, vol. 75, n. 253, 223-229, 2005.

* Jack polynomials. <https://www.symmetricfunctions.com/jack.htm>.
