# fortran-bessels
Fortran port (stub) of the [Bessels.jl](https://github.com/heltonmc/Bessels.jl.git) repository

Building

Currently, only the `besselj0(x)` and `besselj1(x)` functions are available. It is in the `bessel_constants` module. 

A simple build can be achieved by running 

```
gfortran -ffree-line-length-none -O3 -march=native -ffast-math src/bessels_constants.f90 src/bessels.f90 test/bessels_test.f90 -o bessels_test.exe
```

These are the results of a sample performance test on an M1 Mac with gfortran 12.1.0: 

```
[bessel_j0] INTRINSIC time used:   68.2944 ns/eval, sum(z)=49.637464310454263
[bessel_j0] PACKAGE   time used:   32.5004 ns/eval, sum(z)=49.637464310574316
[bessel_j1] INTRINSIC time used:   65.5116 ns/eval, sum(z)=-156.97177970030080
[bessel_j1] PACKAGE   time used:   33.2118 ns/eval, sum(z)=-75.848986598550923
```

this package is approximately *2x faster* than gcc's intrinsic function.


