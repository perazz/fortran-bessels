# Bessel Functions of the First Kind: J {#besselj}

The functions \f$ J_0(x) \f$, \f$ J_1(x) \f$, and \f$ J_n(x) \f$ are the solutions
of Bessel's differential equation

\f[
    x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} + (x^2 - \nu^2)\, y = 0
\f]

that remain finite at the origin. For non-negative integer order \f$ n \f$ they
admit the series

\f[
    J_n(x) = \sum_{m=0}^{\infty} \frac{(-1)^m}{m!\,(m+n)!}
             \left(\frac{x}{2}\right)^{2m+n}.
\f]

\f$ J_\nu(x) \f$ is oscillatory and decays like
\f$ \sqrt{2/(\pi x)}\,\cos\!\left(x - \nu\pi/2 - \pi/4\right) \f$ as
\f$ x \to \infty \f$.

## Function shape

![J0 and J1](besselj.png)

![Integer-order J0 through J4](besseljn.png)

## Domain and special values

| \f$ x \f$        | \f$ J_0(x) \f$         | \f$ J_1(x) \f$          |
|------------------|------------------------|-------------------------|
| \f$ 0 \f$        | \f$ 1 \f$              | \f$ 0 \f$               |
| \f$ -x \f$       | \f$ J_0(x) \f$ (even) | \f$ -J_1(x) \f$ (odd)   |
| \f$ +\infty \f$  | \f$ 0 \f$             | \f$ 0 \f$               |

The J family is defined for all real \f$ x \f$. `besseljn(n, x)` follows the
parity of `n`: even `n` gives an even function, odd `n` an odd function.

## Usage

```fortran
use bessels, only: besselj0, besselj1, besseljn, BK

real(BK) :: x(3), y(3)
x = [0.0_BK, 1.0_BK, 2.5_BK]

y = besselj0(x)               ! elemental over the whole array
print *, besselj1(2.5_BK)
print *, besseljn(3, 2.5_BK)  ! integer order n = 3
```

## Accuracy and performance

`besselj0`/`besselj1` use minimax polynomials near the origin, root- and
extremum-centred polynomials on the middle range, and a Hankel asymptotic
expansion for large \f$ x \f$ (the `SQ2OPI`/`PIO4` phase form). On the benchmark
machine they run roughly **2× faster** than the gfortran intrinsics
`bessel_j0`/`bessel_j1` while matching them to within a few ULP. `besseljn`
composes the fixed-order kernels through a stable recurrence.

The figure below shows the relative error of this library against an `mpmath`
arbitrary-precision reference, in ULPs of `real64`:

![J0 / J1 relative error vs reference](besselj_accuracy.png)

@see @ref bessely
@see @ref besseli
@see @ref besselk
@see @ref hankel
@see @ref mainpage
