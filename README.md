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
 gfortran -ffree-line-length-none -O3 -march=native -ffast-math src/3rd_party/ribesl.f src/3rd_party/rkbesl.f src/bessels_constants.f90 src/bessels.f90 test/bessels_test.f90 -o bessels_test.exe
```

These are the results of a sample performance test on an M1 Mac with gfortran 12.1.0.
For the functions where an intrinsic Fortran equivalent is available, the intrinsic version is compared against.
For all others, the [netlib specfun](https://netlib.org/specfun/) package is employed, in the current refactoring by [Scivision](https://github.com/scivision/rpn-calc-fortran).
For the non-fortran-intrinsic functions, this package is ludicrously faster than the netlib counterpart!

```
[bessel_j0] INTRINSIC time used:   37.5113 ns/eval, sum(z)=9476.3324505667606
[bessel_j0] PACKAGE   time used:   17.8369 ns/eval, sum(z)=9476.3324505666478
[bessel_j1] INTRINSIC time used:   36.7986 ns/eval, sum(z)=-284.46168826127564
[bessel_j1] PACKAGE   time used:   17.8452 ns/eval, sum(z)=-284.46168826129275
[bessel_y0] INTRINSIC time used:   28.4847 ns/eval, sum(z)=1376.4176554633455
[bessel_y0] PACKAGE   time used:   18.0247 ns/eval, sum(z)=1376.4176554633682
[bessel_y1] INTRINSIC time used:   28.5509 ns/eval, sum(z)=-33903.574400809193
[bessel_y1] PACKAGE   time used:   17.9210 ns/eval, sum(z)=-33903.574400809302
[bessel_k0] NETLIB    time used:   44.2205 ns/eval, sum(z)=168876.38538504631
[bessel_k0] PACKAGE   time used:    6.0067 ns/eval, sum(z)=168876.38538504628
[bessel_k1] NETLIB    time used:   27.2245 ns/eval, sum(z)=29117.807091784642
[bessel_k1] PACKAGE   time used:    5.9314 ns/eval, sum(z)=448921.45244578301
[bessel_i0] NETLIB    time used: 1962.9280 ns/eval, sum(z)=0.95961921716826134E+263
[bessel_i0] PACKAGE   time used:   10.3035 ns/eval, sum(z)=0.95961921716826120E+263
[bessel_i1] NETLIB    time used:  479.3809 ns/eval, sum(z)=0.11073899685120145E+48
[bessel_i1] PACKAGE   time used:   10.4763 ns/eval, sum(z)=0.11017821571878319E+48

```

this package is approximately *2x faster* than gcc's intrinsic function.


