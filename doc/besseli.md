# Modified Bessel Functions of the First Kind: I {#besseli}

The functions \f$ I_0(x) \f$ and \f$ I_1(x) \f$ solve the modified Bessel equation

\f[
    x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} - (x^2 + \nu^2)\, y = 0
\f]

and grow exponentially with \f$ x \f$. They are obtained from the ordinary Bessel
functions through \f$ I_\nu(x) = i^{-\nu} J_\nu(i x) \f$, and admit the series

\f[
    I_n(x) = \sum_{m=0}^{\infty} \frac{1}{m!\,(m+n)!}
             \left(\frac{x}{2}\right)^{2m+n},
\f]

with the large-\f$ x \f$ asymptotic
\f$ I_\nu(x) \sim e^{x} / \sqrt{2\pi x} \f$.

## Function shape

![I0 and I1](besseli.png)

## Domain and special values

| \f$ x \f$        | \f$ I_0(x) \f$         | \f$ I_1(x) \f$          |
|------------------|------------------------|-------------------------|
| \f$ 0 \f$        | \f$ 1 \f$              | \f$ 0 \f$               |
| \f$ -x \f$       | \f$ I_0(x) \f$ (even) | \f$ -I_1(x) \f$ (odd)   |
| \f$ +\infty \f$  | \f$ +\infty \f$       | \f$ +\infty \f$         |

Both functions are defined for all real \f$ x \f$ (the implementation works on
\f$ |x| \f$, restoring the sign for the odd \f$ I_1 \f$).

## Usage

```fortran
use bessels, only: besseli0, besseli1, BK

real(BK) :: x(3), y(3)
x = [0.0_BK, 1.0_BK, 3.0_BK]

y = besseli0(x)            ! elemental over the array
print *, besseli1(3.0_BK)
```

## Accuracy and performance

`besseli0`/`besseli1` switch at \f$ x = 7.75 \f$: a minimax polynomial in
\f$ (x/2)^2 \f$ below, and a scaled \f$ \sqrt{x}\,e^{-x} I_\nu(x) \f$ minimax form
(Remez-fitted) above. Compared with the netlib `RIBESL` reference they are
**dramatically faster** (hundreds of ns/eval → ~10 ns/eval) at matching accuracy:

![I0 / I1 relative error vs reference](besseli_accuracy.png)

@see @ref besselk
@see @ref besselj
@see @ref mainpage
