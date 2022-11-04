# fortran-bessels
Fortran port (stub) of the [Bessels.jl](https://github.com/heltonmc/Bessels.jl.git) repository

Building

Currently available functions are in the `bessel_constants` module: 
- `besselj0(x)`
- `besselj1(x)`
- `bessely0(x)`
- `bessely1(x)`

A simple build can be achieved by running: 

```
gfortran -ffree-line-length-none -O3 -march=native -ffast-math src/bessels_constants.f90 src/bessels.f90 test/bessels_test.f90 -o bessels_test.exe
```

These are the results of a sample performance test on an M1 Mac with gfortran 12.1.0: 

```
[bessel_j0] INTRINSIC time used:   37.5113 ns/eval, sum(z)=9476.3324505667606
[bessel_j0] PACKAGE   time used:   17.8369 ns/eval, sum(z)=9476.3324505666478
[bessel_j1] INTRINSIC time used:   36.7986 ns/eval, sum(z)=-284.46168826127564
[bessel_j1] PACKAGE   time used:   17.8452 ns/eval, sum(z)=-284.46168826129275
[bessel_y0] INTRINSIC time used:   28.4847 ns/eval, sum(z)=1376.4176554633455
[bessel_y0] PACKAGE   time used:   18.0247 ns/eval, sum(z)=1376.4176554633682
[bessel_y1] INTRINSIC time used:   28.5509 ns/eval, sum(z)=-33903.574400809193
[bessel_y1] PACKAGE   time used:   17.9210 ns/eval, sum(z)=-33903.574400809302

```

this package is approximately *2x faster* than gcc's intrinsic function.


