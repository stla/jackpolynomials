1.0.0.0
-------
* initial release

1.0.0.1
-------
* removed the upper bounds of the dependencies

1.1.0.0
-------
* replaced the 'mpolynomials' dependency with 'hspray'
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
* it is now possible to choose which Jack polynomial to get or evaluate, `J`, `C`, `P` or `Q` 
(the previous versions returned `J` only)

* it is now possible to get Jack polynomials with a symbolic Jack parameter