# jackpolynomials

***Jack, zonal, Schur, and other symmetric polynomials.***

<!-- badges: start -->
[![Stack-lts](https://github.com/stla/jackpolynomials/actions/workflows/Stack-lts.yml/badge.svg)](https://github.com/stla/jackpolynomials/actions/workflows/Stack-lts.yml)
[![Stack-nightly](https://github.com/stla/jackpolynomials/actions/workflows/Stack-nightly.yml/badge.svg)](https://github.com/stla/jackpolynomials/actions/workflows/Stack-nightly.yml)
<!-- badges: end -->

Schur polynomials have applications in combinatorics and zonal polynomials have
applications in multivariate statistics. They are particular cases of
[Jack polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to compute these polynomials. It also allows to compute other 
symmetric polynomials: Kostka-Foulkes polynomials, t-Schur polynomials, 
Hall-Littlewood polynomials, Kostka-Macdonald polynomials, and Macdonald 
polynomials.

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
import Math.Algebra.SymmetricPolynomials
jp = jackPol' 3 [3, 1, 1] 2 'J'
putStrLn $ prettySymmetricQSpray jp
-- 42*M[3,1,1] + 28*M[2,2,1]
```

And another example, with a symbolic Jack polynomial:

```haskell
import Math.Algebra.JackSymbolicPol
import Math.Algebra.SymmetricPolynomials
jp = jackSymbolicPol' 3 [3, 1, 1] 'J'
putStrLn $ prettySymmetricParametricQSpray ["a"] jp
-- { [ 4*a^2 + 10*a + 6 ] }*M[3,1,1] + { [ 8*a + 12 ] }*M[2,2,1]
```

Of course you can use these functions for other polynomials, but carefully: 
they do not check the symmetry. This new module provides the function 
`isSymmetricSpray` to check the symmetry of a polynomial, much more efficient 
than the function with the same name in the **hspray** package.


### Hall inner product

As of version 1.4.1.0, the package provides an implementation of the Hall 
inner product with parameter. It is known that the Jack polynomials with 
Jack parameter $\alpha$ are orthogonal for the Hall inner product with 
parameter $\alpha$. 

There is a function `hallInnerProduct` as well as a function 
`symbolicHallInnerProduct`. The latter allows to get the Hall inner product 
of two symmetric polynomials without substituting a value to the parameter 
$\alpha$. The Hall inner product of two symmetric polynomials is a polynomial 
in $\alpha$, so the result of `symbolicHallInnerProduct` is a `Spray` object.

Let's see a first example with a power sum polynomial. These symmetric 
polynomials are implemented in the package. We display the result by using 
`alpha` to denote the parameter of the Hall product.

```haskell
import Math.Algebra.SymmetricPolynomials 
import Math.Algebra.Hspray hiding (psPolynomial)
psPoly = psPolynomial 4 [2, 1, 1] :: QSpray
hip = symbolicHallInnerProduct psPoly psPoly
putStrLn $ prettyQSprayXYZ ["alpha"] hip
-- 4*alpha^3
```

Now let's consider the following situation. We want to get the symbolic Hall 
inner product of a Jack polynomial with itself, and we deal with a symbolic 
Jack parameter in this polynomial. We denote it by `t` to distinguish it from 
the parameter of the Hall product that we still denote by `alpha`.

The signature of the `symbolicHallInnerProduct` is a bit misleading:
```haskell
Spray a -> Spray a -> Spray a
```
because the `Spray a` of the output is not of the same family as the two 
`Spray a` inputs: this is a univariate polynomial in $\alpha$. 

We use the function `jackSymbolicPol'` to compute a Jack polynomial. It
returns a `ParametricQSpray` spray, a type alias of `Spray RatioOfQSprays`.

```haskell
import Math.Algebra.JackSymbolicPol
import Math.Algebra.SymmetricPolynomials 
import Math.Algebra.Hspray 
jp = jackSymbolicPol' 2 [3, 1] 'P'
hip = symbolicHallInnerProduct jp jp
putStrLn $ prettyParametricQSprayABCXYZ ["t"] ["alpha"] hip
-- { [ 3*t^2 + 6*t + 11 ] %//% [ t^2 + 2*t + 1 ] }*alpha^2 + { [ 4*t^2 + 16*t + 16 ] %//% [ t^2 + 2*t + 1 ] }*alpha
```

One could be interested in computing the Hall inner product of a Jack 
polynomial with itself when the Jack parameter and the parameter of the 
Hall product are identical. That is, we want to take `alpha = t` in the 
above expression. Since the symbolic Hall product is a `ParametricQSpray` 
spray, one can substitute its variable `alpha` by a `RatioOfQSprays` 
object. On the other hand, `t` represents a `QSpray` object, but one can 
identify a `QSpray` to a `RatioOfQSprays` by taking the unit spray as the 
denominator, that is, by applying the `asRatioOfSprays` function. Finally 
we get the desired result if we evaluate the symbolic Hall product by 
replacing `alpha` with `asRatioOfSprays (qlone 1)`, since `t` is the 
first polynomial variable, `qlone 1`. 

```haskell
prettyRatioOfQSpraysXYZ ["t"] $ evaluate hip [asRatioOfSprays (qlone 1)]
-- [ 3*t^4 + 10*t^3 + 27*t^2 + 16*t ] %//% [ t^2 + 2*t + 1 ]
```


### Hall-Littlewood polynomials

The package can also compute the Hall-Littlewood polynomials. A Hall-Littlewood 
polynomial is a multivariate symmetric polynomial associated to an integer 
partition and whose coefficients depend on a parameter. More precisely, the 
coefficients are some polynomials in this parameter. So the Hall-Littlewood 
polynomials implemented in the package, returned by the 
`hallLittlewoodPolynomial` function, are represented by some sprays of type 
`SimpleParametricSpray a`, an alias of the type `Spray (Spray a)`. 

When the value of the parameter of a Hall-Littlewood polynomial is `0`, then 
this polynomial is the Schur polynomial of the given partition.

```haskell
import Math.Algebra.JackPol
import Math.Algebra.SymmetricPolynomials 
import Math.Algebra.Hspray 
lambda = [2, 1]
hlPoly = hallLittlewoodPolynomial' 3 lambda 'P' 
putStrLn $ prettySymmetricSimpleParametricQSpray ["t"] hlPoly
-- (1)*M[2,1] + (-t^2 - t + 2)*M[1,1,1]
hlPolyAt0 = substituteParameters hlPoly [0]
hlPolyAt0 == schurPol' 3 lambda
-- True
```


### Macdonald polynomials

As of version 1.4.5.0, the package can compute some Macdonald polynomials. 
The Macdonald polynomials are symmetric multivariate polynomials whose 
coefficients depend on two parameters usually denoted by $q$ and $t$. 
Let's consider for example the Macdonald P-polynomial. It generalizes 
the Hall-Littlewood P-polynomial: the Hall-Littlewood P-polynomial with 
parameter $t$ is obtained from the Macdonald P-polynomial by substituting
the parameter $q$ with $0$.

```haskell
import Math.Algebra.Hspray
import Math.Algebra.SymmetricPolynomials
n = 4
lambda = [2, 2]
macPoly = macdonaldPolynomial' n lambda 'P'
poly = changeParameters macPoly [zeroSpray, qlone 1]
hlPoly = hallLittlewoodPolynomial' n lambda 'P'
asSimpleParametricSpray poly == hlPoly
-- True
```

We use `asSimpleParametricSpray` because, contrary to the Hall-Littlewood 
J-polynomial, the Macdonald J-polynomial is not represented by a 
`SimpleParametricSpray a` spray but by a `ParametricSpray a` spray, because 
its coefficients are not polynomials in the two parameters $q$ and $t$, but 
ratios of polynomials.


## References

* I.G. Macdonald. *Symmetric Functions and Hall Polynomials*. Oxford Mathematical Monographs. The Clarendon Press Oxford University Press, New York, second edition, 1995.

* J. Demmel and P. Koev. *Accurate and efficient evaluation of Schur and Jack functions*. Mathematics of computations, vol. 75, n. 253, 223-229, 2005.

* Jack polynomials. <https://www.symmetricfunctions.com/jack.htm>.
