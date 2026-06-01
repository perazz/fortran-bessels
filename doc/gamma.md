# Gamma Function {#gamma}

`gamma_BK(x)` evaluates Euler's gamma function

\f[
    \Gamma(x) = \int_0^{\infty} t^{x-1} e^{-t}\, dt,
\f]

the continuous extension of the factorial, \f$ \Gamma(n) = (n-1)! \f$ for positive
integers \f$ n \f$. It is used internally by the variable-order Bessel routines and
is re-exported from the `bessels` module for convenience.

## Function shape

![Gamma function](gamma.png)

## Domain and special values

| \f$ x \f$                  | \f$ \Gamma(x) \f$                 |
|----------------------------|-----------------------------------|
| \f$ 1 \f$                  | \f$ 1 \f$                         |
| positive integer \f$ n \f$ | \f$ (n-1)! \f$                    |
| \f$ 0, -1, -2, \dots \f$   | simple poles (\f$ \pm\infty \f$) |
| \f$ +\infty \f$            | \f$ +\infty \f$                   |

`gamma_BK` accepts a real argument; an integer-argument overload provides a
factorial-optimized fast path for small integer orders.

## Usage

```fortran
use bessels, only: gamma_BK, BK

real(BK) :: x(3), g(3)
x = [0.5_BK, 1.0_BK, 5.0_BK]

g = gamma_BK(x)          ! elemental:  [sqrt(pi), 1, 24]
```

## Accuracy and performance

`gamma_BK` is faster than the gfortran intrinsic `gamma` (≈26 ns/eval vs
≈38 ns/eval) at matching accuracy. The relative error against an `mpmath`
arbitrary-precision reference is shown below in ULPs of `real64`:

![Gamma relative error vs reference](gamma_accuracy.png)

@see @ref besselj
@see @ref mainpage
