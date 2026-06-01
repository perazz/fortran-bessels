# Hankel Functions {#hankel}

The Hankel functions (Bessel functions of the third kind) are the complex
combinations

\f[
    H^{(1)}_\nu(x) = J_\nu(x) + i\,Y_\nu(x), \qquad
    H^{(2)}_\nu(x) = J_\nu(x) - i\,Y_\nu(x),
\f]

which behave like outgoing and incoming cylindrical waves,
\f$ H^{(1,2)}_\nu(x) \sim \sqrt{2/(\pi x)}\,
\exp\!\left[\pm i\left(x - \nu\pi/2 - \pi/4\right)\right] \f$ as
\f$ x \to \infty \f$. For real \f$ x \f$ and real \f$ \nu \f$ they are complex
conjugates, \f$ H^{(2)}_\nu(x) = \overline{H^{(1)}_\nu(x)} \f$.

## Function shape

![Real and imaginary parts of the first-kind Hankel function](hankel.png)

## Domain and special values

| Quantity | Value |
|----------|-------|
| `besselh(nu, k, x)` | \f$ H^{(k)}_\nu(x) \f$, returns `complex(BK)`, with `k` ∈ {1, 2} |
| `hankelh1(nu, x)`   | \f$ H^{(1)}_\nu(x) \f$ = `besselh(nu, 1, x)`                     |
| `hankelh2(nu, x)`   | \f$ H^{(2)}_\nu(x) \f$ = `besselh(nu, 2, x)`                     |
| \f$ x \le 0 \f$     | `NaN + NaN·i` (built on \f$ Y_\nu \f$, undefined there)         |

## Usage

```fortran
use bessels, only: besselh, hankelh1, hankelh2, BK

complex(BK) :: h

h = hankelh1(0.0_BK, 5.0_BK)         ! H^(1)_0(5)
print *, h, hankelh2(0.0_BK, 5.0_BK) ! conjugate for real argument
print *, besselh(1.0_BK, 2, 5.0_BK)  ! H^(2)_1(5)
```

@warning `besselh`, `hankelh1`, and `hankelh2` are validated for **integer**
order \f$ \nu \f$ only. They are evaluated by composing \f$ J_\nu + iY_\nu \f$,
so they inherit the non-integer-\f$ \nu \f$ bugs of @ref bessely. Non-integer
\f$ \nu \f$ will return incorrect values; see
[roadmap item 07](https://github.com/perazz/fortran-bessels/blob/main/todo/07-nonintegerNU-bugs.md).

## Implementation

`besselh` composes the @ref besselj and @ref bessely kernels directly. The fast
large-\f$ x \f$ `hankel_debye` path is currently bypassed pending the
non-integer-\f$ \nu \f$ fix, so accuracy follows that of \f$ J_\nu \f$ and
\f$ Y_\nu \f$ for integer orders.

@see @ref besselj
@see @ref bessely
@see @ref mainpage
