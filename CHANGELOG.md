1.0.0.0
-------
* initial release

1.0.0.1
-------
* removed the upper bounds of the dependencies

1.1.0.0
-------
* replaced the **mpolynomials** dependency with **hspray**
* unit tests

1.1.0.1
-------
* unexported some useless functions
* one more unit test

1.1.1.0
-------
* `schurPol` now returns a `Spray a`
* added package upper bounds in the cabal file
* increased the version of the dependencies **hspray** and **hypergeomatrix**
* cleaned the code
* tested with higher versions of GHC
* new unit tests

1.1.2.0
-------
* skew Schur polynomials (functions `skewSchur` and `skewSchurPol`)

1.2.0.0
-------
* it is now possible to choose which Jack polynomial to get or evaluate, 
`J`, `C`, `P` or `Q` (the previous versions returned `J` only)

* it is now possible to get Jack polynomials with a symbolic Jack parameter

1.2.1.0
-------
* a new module provides some stuff to deal with symmetric polynomials, mainly 
some functions to print them as a linear combination of the monomial symmetric 
polynomials, and a function to check the symmetry

1.2.2.0
-------
* slight modifications due to the upgrade of **hspray**

1.3.0.0
-------
* the type of the Jack polynomials with a symbolic Jack parameter has changed 
from `OneParameterSpray a` to `ParametricSpray a`

1.4.0.0
-------
* upgraded the **hspray** dependency (an error has been fixed in this new 
version)

* added the Laplace-Beltrami operator and the Calogero-Sutherland operator; 
the Jack polynomials are eigenpolynomials of these operators

1.4.1.0
-------
* new function `psPolynomial`, to get a power sum symmetric polynomial

* new function `psCombination`, to get a symmetric polynomial as a linear 
combination of some power sum polynomials

* new function `hallInnerProduct`, to compute the Hall inner product between 
two symmetric polynomials, aka the Jack-scalar product or the deformed Hall 
inner product; there is also the function `symbolicHallInnerProduct`, to get 
the Hall inner product with a symbolic parameter

1.4.2.0
-------
* new function `cshPolynomial`, to get a complete symmetric homogeneous polynomial

* new function `cshCombination`, to get a symmetric polynomial as a linear 
combination of some complete symmetric homogeneous polynomials

* new function `esPolynomial`, to get an elementary symmetric polynomial

* new function `esCombination`, to get a symmetric polynomial as a linear 
combination of some elementary symmetric polynomials

* new function `schurCombination`, to get a symmetric polynomial as a linear 
combination of some Schur polynomials

* new function `jackCombination`, to get a symmetric polynomial as a linear 
combination of some Jack polynomials with a fixed Jack parameter

* new function `jackSymbolicCombination`, to get a symmetric polynomial as a linear 
combination of some Jack polynomials with symbolic Jack parameter

* new functions `kostkaNumbers` and `symbolicKostkaNumbers`, to get the Kostka 
numbers with parameter

1.4.3.0
-------
* new function `kostkaFoulkesPolynomial`, to get a Kostka-Foulkes polynomial

* new function `hallLittlewoodPolynomial`, to get a Hall-Littlewood polynomial

* new function `skewHallLittlewoodPolynomial`, to get a skew Hall-Littlewood 
polynomial

* new function `flaggedSchurPol`, to get a flagged Schur polynomial