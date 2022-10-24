# fortran-bessels
Fortran port (stub) of the [Bessels.jl](https://github.com/heltonmc/Bessels.jl.git) repository

Building

Currently, only the `besselj0(x)` function is available. It is in the `bessel_constants` module. 

A simple build can be achieved by running 

```
gfortran -O3 -march=native -ffast-math src/bessels_constants.f90 src/bessels.f90 test/bessels_test.f90 -o bessels_test.exe
```



