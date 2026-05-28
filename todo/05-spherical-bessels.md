# 05 — Spherical Bessel functions

## Goal

Public `sphericalbesselj(nu, x)`, `sphericalbessely(nu, x)`, `sphericalbesseli(nu, x)`, `sphericalbesselk(nu, x)` — all `elemental real(BK)`.

These are half-integer-order rescalings of the regular Bessels:
- `j_n(x) = √(π/(2x)) · J_{n+1/2}(x)`
- `y_n(x) = √(π/(2x)) · Y_{n+1/2}(x)`
- `i_n(x) = √(π/(2x)) · I_{n+1/2}(x)`
- `k_n(x) = √(2/(πx)) · K_{n+1/2}(x)`  (note: different prefactor)

For integer `n` they have closed forms in elementary functions for small `n` (which is most of the practical use).

## Dependencies

- `sphericalbesselj/y` need `besselj(nu, x)` / `bessely(nu, x)` for non-integer `n` — see [01](01-wire-existing-variable-order.md).
- `sphericalbesseli` needs `besseli(nu, x)` — [04](04-besseli-variable-order.md).
- `sphericalbesselk` needs `besselk(nu, x)` — [03](03-besselk-variable-order.md).

For **integer** `n` we can sidestep the half-integer Bessel dependency entirely via the closed forms below — so if [03] and [04] are delayed, an integer-only version of all four spherical functions ships first.

## What's already there

Nothing spherical. All the dependencies are tracked in other tier files.

## Approach — port [Bessels.jl/src/BesselFunctions/sphericalbessel.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/sphericalbessel.jl) and [modifiedsphericalbessel.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/modifiedsphericalbessel.jl)

### Integer `n` fast paths (no Bessel dependency)

Base cases:
- `j_0(x) = sin(x)/x`, `j_1(x) = sin(x)/x² − cos(x)/x`
- `y_0(x) = −cos(x)/x`, `y_1(x) = −cos(x)/x² − sin(x)/x`
- `i_0(x) = sinh(x)/x`, `i_1(x) = (x cosh(x) − sinh(x))/x²`
- `k_0(x) = e^{−x}/x`, `k_1(x) = (1/x + 1/x²) · e^{−x}`

Recurrence (same shape for all four):
`f_{n+1}(x) = (2n+1)/x · f_n(x) − f_{n−1}(x)`  (for j, y, modified-i)
`f_{n+1}(x) = (2n+1)/x · f_n(x) + f_{n−1}(x)`  (for k — sign flip)

**Stability:** forward recurrence is stable for `y_n` (any `x`) and for `k_n` (any `x`), and for `j_n` / `i_n` when `n < x` — same considerations as the regular Bessels. For `j_n` with `n > x` use either backward Miller recurrence or fall through to `besselj(n+1/2, x)` once [01](01-wire-existing-variable-order.md) is done.

### Cutoffs per Bessels.jl

- Integer `n < 250` and `x ≥ n`: forward recurrence.
- Integer `n < 60` and `x < n`: combine forward `y_n` recurrence with a continued-fraction `j_n/j_{n+1}` and the Wronskian `j_n y_{n+1} − j_{n+1} y_n = 1/x²` to recover `j_n`.
- Otherwise: half-integer reduction via `besselj_positive_args(n + 0.5, x)`.
- Very small `x`: power series `j_n(x) = x^n/(2n+1)!! · (1 − x²/(2(2n+3)) + …)` to avoid `sin(x)/x` cancellation.

### Modified spherical (i, k) — [modifiedsphericalbessel.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/modifiedsphericalbessel.jl)

- `i_n`: closed-form `sinh/cosh` expressions for `n ∈ {0,1,2}`, power series for small `x`, else half-integer reduction through `besseli(n+1/2, x)`.
- `k_n`: forward recurrence from explicit `k_0`, `k_1` for `n < 41.5`. Otherwise `√(π/(2x)) · K_{n+1/2}(x)`. Negative `n` maps via `k_{−n}(x) = k_{n−1}(x)`.

### Real `nu` (non-integer)

Always reduces to the half-integer rescaling:
```
sphericalbesselj(nu, x) = sqrt(PIO2/x) * besselj(nu + HALF, x)
```
…and analogously for the other three. No new code beyond the prefactor.

## Step-by-step

1. **New module** [src/bessels_spherical.f90](../src/bessels_spherical.f90), `use bessels_constants`, `use bessels` (for `besselj_positive_args`, `bessely_positive_args`, `besseli`, `besselk` once they're public — make them so).
2. **`sphericalbesselj_int(n, x)`** — closed-form `n=0`, recurrence for `n < min(250, x)`, continued-fraction + Wronskian for `n > x` and `n < 60`, else half-integer reduction.
3. **`sphericalbessely_int(n, x)`** — closed-form `n=0,1`, forward recurrence from `y_0, y_1`.
4. **`sphericalbesseli_int(n, x)`, `sphericalbesselk_int(n, x)`** — same shape.
5. **Public elemental `sphericalbesselj(nu, x)`** dispatches: if `nu == real(int(nu), BK) .and. nu >= 0` → integer fast path; else `sqrt(PIO2/x) * besselj_positive_args(nu + HALF, x)` (with sign/domain handling for `x < 0` via parity).
6. **Domain.** `x = 0` limits: `j_0(0) = 1`, `j_{n≥1}(0) = 0`, `y_n(0) = -∞` for all `n`, `i_0(0) = 1`, `i_{n≥1}(0) = 0`, `k_n(0) = +∞`. Match the Bessels.jl values exactly.
7. **Wire into `bessels.f90`.** Export the four functions.

## Tests

- `test_sphericalbesselj` — `n ∈ {0, 1, 5, 50}`, `x ∈ {1.0, 10.0, 100.0}` against tabulated reference and against `sin(x)/x`-style closed forms for low `n`.
- `test_sphericalbessely_recurrence_stability` — `n ∈ {0..60}`, `x = 1.0` — verify `y_n` grows monotonically and never NaN.
- `test_sphericalbesseli`, `test_sphericalbesselk` — analogous.
- `test_sphericalbessel_half_integer_consistency` — verify `sphericalbesselj(5.5, x) == sqrt(PIO2/x) * besselj(6.0, x)` to `eps(BK) * 10` (this checks the dispatcher routes through `besselj_positive_args` correctly).
- Benchmarks for all four.

## Done when

- Four spherical functions exported.
- Integer-`n` fast paths verified independent of the dependencies (they don't require `besselj/y/i/k` of fractional order).
- Half-integer reduction tested once those dependencies are public.
- README "Not yet implemented" loses four entries.
