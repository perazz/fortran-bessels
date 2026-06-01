# Modified Bessel Functions of the Second Kind: K {#besselk}

The functions \f$ K_0(x) \f$ and \f$ K_1(x) \f$ are the second solutions of the
modified Bessel equation

\f[
    x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} - (x^2 + \nu^2)\, y = 0,
\f]

decaying exponentially with \f$ x \f$. They can be written as

\f[
    K_\nu(x) = \frac{\pi}{2}\,\frac{I_{-\nu}(x) - I_\nu(x)}{\sin(\nu\pi)},
\f]

with the large-\f$ x \f$ asymptotic
\f$ K_\nu(x) \sim \sqrt{\pi/(2x)}\, e^{-x} \f$.

## Function shape

![K0 and K1](besselk.png)

## Domain and special values

| \f$ x \f$          | \f$ K_0(x) \f$ / \f$ K_1(x) \f$            |
|--------------------|--------------------------------------------|
| \f$ x \le 0 \f$    | `NaN` (undefined for \f$ x<0 \f$, diverges as \f$ x\to 0^+ \f$) |
| \f$ x \to 0^+ \f$  | \f$ +\infty \f$                            |
| \f$ +\infty \f$    | \f$ 0 \f$                                  |

## Usage

```fortran
use bessels, only: besselk0, besselk1, BK

real(BK) :: x(2), y(2)
x = [0.5_BK, 2.0_BK]

y = besselk0(x)            ! elemental over the array
print *, besselk1(2.0_BK)
```

## Accuracy and performance

`besselk0`/`besselk1` use two branches with rational (boost-style, Holoborodko)
approximations: \f$ x \le 1 \f$ uses the \f$ \log x \cdot I_\nu(x) \f$ correction
form, and \f$ x > 1 \f$ a scaled \f$ \sqrt{x}\,e^{x} K_\nu(x) \f$ rational form.
Against the netlib `RKBESL` reference they are several times faster
(~6 ns/eval vs ~27–44 ns/eval) at matching accuracy:

![K0 / K1 relative error vs reference](besselk_accuracy.png)

@see @ref besseli
@see @ref bessely
@see @ref mainpage
