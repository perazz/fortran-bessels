# Bessel Functions of the Second Kind: Y {#bessely}

The functions \f$ Y_0(x) \f$, \f$ Y_1(x) \f$, and \f$ Y_\nu(x) \f$ are the second,
linearly independent solutions of Bessel's equation

\f[
    x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} + (x^2 - \nu^2)\, y = 0,
\f]

singular at the origin. They are related to the first kind by

\f[
    Y_\nu(x) = \frac{J_\nu(x)\cos(\nu\pi) - J_{-\nu}(x)}{\sin(\nu\pi)},
\f]

and share the large-\f$ x \f$ asymptotic
\f$ Y_\nu(x) \sim \sqrt{2/(\pi x)}\,\sin\!\left(x - \nu\pi/2 - \pi/4\right) \f$.

## Function shape

![Y0 and Y1](bessely.png)

## Domain and special values

| \f$ x \f$          | \f$ Y_0(x) \f$ / \f$ Y_1(x) \f$           |
|--------------------|-------------------------------------------|
| \f$ x < 0 \f$      | `NaN` (undefined for real \f$ x<0 \f$)    |
| \f$ x \to 0^+ \f$  | \f$ -\infty \f$ (returned as `-huge(BK)`) |
| \f$ +\infty \f$    | \f$ 0 \f$                                 |

## Usage

```fortran
use bessels, only: bessely0, bessely1, bessely, BK

real(BK) :: x(2), y(2)
x = [1.0_BK, 7.5_BK]

y = bessely0(x)               ! elemental over the array
print *, bessely1(7.5_BK)
print *, bessely(2.0_BK, 7.5_BK)   ! real-order Yν (integer ν here)
```

@warning `bessely(nu, x)` is validated for **integer** order \f$ \nu \f$ only.
Non-integer \f$ \nu \f$ currently routes through variable-order paths with known
port bugs (power-series argument swap, Chebyshev mapping, `hankel_debye` complex
output) and will return incorrect values. See
[roadmap item 07](https://github.com/perazz/fortran-bessels/blob/main/todo/07-nonintegerNU-bugs.md).

## Accuracy and performance

`bessely0`/`bessely1` use rational approximations on \f$ x \le 5 \f$ (with the
\f$ (2/\pi)\log x \, J_\nu(x) \f$ correction term), Hankel phase/amplitude forms
on the middle range, and a sine asymptotic for large \f$ x \f$. They run roughly
**2× faster** than the gfortran intrinsics `bessel_y0`/`bessel_y1` while matching
them to a few ULP. The accuracy plot below is scoped to \f$ Y_0 \f$ / \f$ Y_1 \f$:

![Y0 / Y1 relative error vs reference](bessely_accuracy.png)

@see @ref besselj
@see @ref hankel
@see @ref besselk
@see @ref mainpage
