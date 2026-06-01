# fortran-bessels {#mainpage}

**fortran-bessels** is a Modern Fortran (2008+) port of
[Bessels.jl](https://github.com/heltonmc/Bessels.jl), built for one thing:
evaluating Bessel functions **as fast and as accurately as possible** in pure
Fortran. The fixed-order kernels beat the gfortran intrinsics by roughly 2× and
the netlib `specfun` reference by orders of magnitude, while matching them to a
few units in the last place (ULP).

Every function is an `elemental real(BK) function` (`BK = real64`), so it accepts
a scalar or a whole array interchangeably — no wrapper loops, no temporaries.

## Function reference

### First kind
- @ref besselj &mdash; \f$ J_0,\; J_1,\; J_n \f$ — oscillatory, finite at the origin

### Second kind
- @ref bessely &mdash; \f$ Y_0,\; Y_1,\; Y_\nu \f$ — oscillatory, singular at the origin

### Modified, first kind
- @ref besseli &mdash; \f$ I_0,\; I_1 \f$ — exponentially growing

### Modified, second kind
- @ref besselk &mdash; \f$ K_0,\; K_1 \f$ — exponentially decaying

### Hankel
- @ref hankel &mdash; \f$ H^{(1)}_\nu,\; H^{(2)}_\nu \f$ — complex-valued

### Support
- @ref gamma &mdash; \f$ \Gamma(x) \f$

## Quick start

```fortran
use bessels, only: besselj0, besselk1, BK

real(BK) :: x(4)
x = [0.5_BK, 1.0_BK, 2.0_BK, 4.0_BK]

print *, besselj0(x)        ! elemental: acts on the whole array
print *, besselk1(2.0_BK)
```

## Building

The canonical build is [fpm](https://fpm.fortran-lang.org):

```bash
fpm test --profile release --flag "-march=native"
```

See the project
[README](https://github.com/perazz/fortran-bessels#readme) for the full benchmark
table and the manual `gfortran` recipe.

@warning The variable-order routines `bessely(nu,x)`, `besselh`, `hankelh1`, and
`hankelh2` are validated for **integer** order \f$ \nu \f$ only. Non-integer order
has known port bugs in the variable-order Y/Hankel paths and will return incorrect
values. See
[roadmap item 07](https://github.com/perazz/fortran-bessels/blob/main/todo/07-nonintegerNU-bugs.md).

## Roadmap

Still to come: real-order \f$ I_\nu,\; K_\nu \f$, spherical Bessels, and the Airy
functions. See the
[implementation roadmap](https://github.com/perazz/fortran-bessels/blob/main/todo/README.md).
