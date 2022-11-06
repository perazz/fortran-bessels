# fortran-bessels
Fortran port (stub) of the [Bessels.jl](https://github.com/heltonmc/Bessels.jl.git) repository

Building

Currently available functions are in the `bessel_constants` module: 
- `besselj0(x)`
- `besselj1(x)`
- `bessely0(x)`
- `bessely1(x)`
- `besselk0(x)`
- `besselk1(x)`

A simple build can be achieved by running: 

```
 gfortran -ffree-line-length-none -O3 -march=native -ffast-math src/legacy/rkbesl.f src/bessels_constants.f90 src/bessels.f90 test/bessels_test.f90 -o bessels_test.exe
```

These are the results of a sample performance test on an M1 Mac with gfortran 12.1.0. 
For the functions where an intrinsic Fortran equivalent is available, the intrinsic version is compared against. For all others, the [netlib specfun](https://netlib.org/specfun/) package is employed, in the current refactoring by [Scivision](https://github.com/scivision/rpn-calc-fortran).   

```
[bessel_j0] INTRINSIC time used:   37.5113 ns/eval, sum(z)=9476.3324505667606
[bessel_j0] PACKAGE   time used:   17.8369 ns/eval, sum(z)=9476.3324505666478
[bessel_j1] INTRINSIC time used:   36.7986 ns/eval, sum(z)=-284.46168826127564
[bessel_j1] PACKAGE   time used:   17.8452 ns/eval, sum(z)=-284.46168826129275
[bessel_y0] INTRINSIC time used:   28.4847 ns/eval, sum(z)=1376.4176554633455
[bessel_y0] PACKAGE   time used:   18.0247 ns/eval, sum(z)=1376.4176554633682
[bessel_y1] INTRINSIC time used:   28.5509 ns/eval, sum(z)=-33903.574400809193
[bessel_y1] PACKAGE   time used:   17.9210 ns/eval, sum(z)=-33903.574400809302
[bessel_k0] NETLIB    time used:   24.1788 ns/eval, sum(z)=1479989096.7517223
[bessel_k0] PACKAGE   time used:    5.9340 ns/eval, sum(z)=16538.641694398026
```

this package is approximately *2x faster* than gcc's intrinsic function.


