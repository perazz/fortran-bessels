# fortran-bessels
A Modern Fortran port of the [Bessels.jl](https://github.com/heltonmc/Bessels.jl.git) library for
fast, accurate Bessel function evaluation in pure Fortran.

📖 **API documentation:** https://perazz.github.io/fortran-bessels

## Available functions

All functions live in the `bessels` module and are `elemental real(BK) function`s
(`BK = real64`) unless noted, so they accept scalars or whole arrays interchangeably.

**Fixed order** — hand-tuned minimax/asymptotic branches, the fast path:
- `besselj0(x)`, `besselj1(x)` — J₀, J₁
- `bessely0(x)`, `bessely1(x)` — Y₀, Y₁
- `besseli0(x)`, `besseli1(x)` — I₀, I₁
- `besselk0(x)`, `besselk1(x)` — K₀, K₁

**Integer order:**
- `besseljn(nu, x)` — Jₙ, `nu` integer

**Real order** (`nu`, `x` both real):
- `bessely(nu, x)` — Yᵥ
- `besselh(nu, k, x) -> complex(BK)` — Hankel Hᵥ⁽ᵏ⁾, `k` ∈ {1, 2}
- `hankelh1(nu, x)`, `hankelh2(nu, x) -> complex(BK)` — Hᵥ⁽¹⁾, Hᵥ⁽²⁾

**Support:**
- `gamma_BK(x)` — gamma function (re-exported from `bessels`)
- `cbrt(x)` — cube-root helper
- constants `BK`, `BSIZE`, `ZERO`, `ONE`, `THIRD`

> ⚠️ **Integer ν only for `bessely`/`besselh`/`hankelh1`/`hankelh2`.**
> These are validated for **integer** orders. Non-integer ν has known port bugs in the
> variable-order Y/Hankel paths (power-series argument swap, Chebyshev mapping, `hankel_debye`
> complex output) and will return incorrect values. See
> [todo/07-nonintegerNU-bugs.md](todo/07-nonintegerNU-bugs.md).

### Not yet implemented
- `besseli(nu, x)`, `besselk(nu, x)` — real-order I, K
- `sphericalbesselj/y/i/k(nu, x)` — spherical Bessels
- `airyai(x)`, `airyaiprime(x)`, `airybi(x)`, `airybiprime(x)` — Airy functions

See the [implementation roadmap](todo/README.md) for the priority order.

## Building

The canonical build is [fpm](https://fpm.fortran-lang.org):

```bash
fpm test --profile release --flag "-march=native"
```

Manual gfortran build (matches the benchmark recipe):

```bash
gfortran -ffree-line-length-none -O3 -march=native -ffast-math \
    src/bessels_constants.f90 src/bessels_gamma.f90 src/bessels_debye.f90 \
    src/bessels.f90 \
    test/3rd_party/ribesl.f90 test/3rd_party/rkbesl.f90 \
    test/bessels_test.f90 -o bessels_test
./bessels_test
```

## Performance

These are the results of a sample performance test on an M1 Mac with gfortran 12.1.0.
The table covers the fixed-order kernels and gamma, where an intrinsic Fortran or netlib
reference exists; the variable-order routines (`besseljn`, `bessely`, Hankel) compose these
kernels and inherit their speed. For functions with an intrinsic Fortran equivalent, the
intrinsic version is compared against. For all others, the
[netlib specfun](https://netlib.org/specfun/) package is employed, in the current refactoring by
[Scivision](https://github.com/scivision/rpn-calc-fortran).

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
[gamma]     INTRINSIC time used:   37.9529 ns/eval, sum(z)=0.14440233357737784E+68
[gamma]     PACKAGE   time used:   26.1333 ns/eval, sum(z)=0.14440233357737787E+68

```

this package is approximately *2x faster* than gcc's intrinsic function. For the
non-fortran-intrinsic functions, this package is ludicrously faster than the netlib counterpart!
