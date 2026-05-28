# 04 — `besseli(nu, x)` variable order

## Goal

Public `besseli(nu, x)` and scaled `besselix(nu, x) = exp(-x) · I_ν(x)`, both `elemental real(BK)`.

## Dependency

**Requires [03-besselk-variable-order.md](03-besselk-variable-order.md) to be done first** — `I_ν(x)` for negative `ν` uses the reflection `I_{-ν}(x) = I_ν(x) + (2/π) sin(πν) K_ν(x)`. Without `besselk(nu, x)` this branch is unimplementable.

The K-port also produces a U-polynomial table for Debye that `besseli_debye` reuses — coordinate the file layout so both modules share it.

## What's already there

- `besseli0(x)`, `besseli1(x)` in [src/bessels.f90](../src/bessels.f90).
- `gamma_BK` for the power series.
- `evalpoly`, `muladd`.
- After [03](03-besselk-variable-order.md) lands: `besselk(nu, x)` for the reflection, plus shared U-polynomial table.

## What's missing

1. I-specific Debye uniform asymptotic — same U-polynomial table as `besselk` but with a different prefactor and an additive sign (NIST 10.41.E3 vs 10.41.E4).
2. Power series for small `x` — `I_ν(x) = (x/2)^ν · Σ (x/2)^{2k} / (k! · Γ(ν+k+1))`.
3. Large-x asymptotic — `I_ν(x) ≈ e^x/√(2πx) · P(1/x)`, where the series coefficients are signed versions of K's (related by `(−1)^k`). Can largely be derived from the K port.
4. Public `besseli(nu, x)` dispatcher.
5. `besselix(nu, x)` — scaled variant, the practically important one for large x.

## Approach — port [Bessels.jl/src/BesselFunctions/besseli.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/besseli.jl)

### Negative arguments and orders

Up front in the dispatcher (before any computation):

```fortran
if (nu == real(int(nu), BK) .and. x < ZERO) then
    ! integer nu, negative x: I_n(-x) = (-1)^n I_n(x)
    besseli = besseli(nu, -x) * (1 - 2*mod(int(abs(nu)), 2))
elseif (nu /= real(int(nu), BK) .and. x < ZERO) then
    besseli = ieee_value(besseli, ieee_quiet_nan)   ! domain error
elseif (nu < ZERO) then
    ! I_{-ν}(x) = I_ν(x) + (2/π) sin(πν) K_ν(x)
    besseli = besseli_positive_args(-nu, x) + TWOOPI * sin(-PI*nu) * besselk(-nu, x)
else
    besseli = besseli_positive_args(nu, x)
end if
```

### Positive-`nu` branches

| Branch | Method |
|---|---|
| Large `x` and not too large `nu` | Asymptotic series in `1/x`, leading factor `e^x / √(2πx)`. |
| Large `nu` or moderate `x` with `nu ≳ x` | Debye uniform asymptotic via I-flavored U-polynomial expansion. |
| Default (small `x`, moderate `nu`) | Power series. |
| Integer `nu`, moderate `x` | Forward recurrence from `besseli0`/`besseli1` using `I_{n+1}(x) = I_{n-1}(x) - (2n/x) I_n(x)` (note: this is the *downward* relation; upward is unstable for I, so for moderate `nu < ~x` use forward, otherwise use the asymptotic branches). Bessels.jl uses a hybrid approach — check the Julia for the exact cutoff. |

### Scaled variant `besselix(nu, x)`

In each branch, factor out `e^x` or `e^{-x}` to keep the result `O(1)` for large `x`. The asymptotic branch naturally has `e^x` separable; the Debye branch produces a `nu·η(p)` exponential that combines cleanly with `e^{-x}` for the scaled form. The power series branch should reconstruct `e^{-x}·I_ν(x)` for `besselix` rather than computing `besseli` and multiplying — cheap with `exp(-x)` applied once after the sum.

## Step-by-step

1. **Pre-flight.** Verify [03](03-besselk-variable-order.md) is done and the shared U-polynomial table is in [src/bessels_debye.f90](../src/bessels_debye.f90) or a shared module.
2. **New module** [src/bessels_besseli.f90](../src/bessels_besseli.f90) — parallel to `bessels_besselk`.
3. **Power series.** Port `besseli_power_series`. Cancellation is minimal because all terms have the same sign — Horner-style accumulation is fine.
4. **Large-x asymptotic.** Port `besseli_large_argument`. The coefficient series is the *unsigned* version of K's (`(μ - 1²)(μ - 3²)…` with `+` between consecutive products instead of `-`). Worth a helper `compute_I_K_asymptotic_coefs(nu)` shared with the K module.
5. **Debye expansion.** Port `besseli_debye`. Same U-table as K, different prefactor `exp(ν·η(p)) · √(p / (2π·ν))` where `p = 1/√(1 + (x/ν)²)`, `η(p) = √(1+(x/ν)²) + log((x/ν)/(1 + √(1+(x/ν)²)))`. For the scaled variant, the `exp(ν·η - x)` combines.
6. **Recurrence.** Forward `besseli_up_recurrence` from `besseli0`, `besseli1` for integer `nu`, gated by a `nu < x_cutoff` check (forward recurrence is stable only when `nu ≪ x`; otherwise use downward Miller recurrence or the Debye branch).
7. **Dispatcher.** Wire branches per the table above. Add cutoff helpers `besseli_debye_cutoff`, `besseli_large_arg_cutoff`.
8. **Wire into `bessels.f90`.** `use bessels_besseli`, export `besseli`, `besselix`.

## Tests

- `test_besseli_nu` — `(nu, x)` grid analogous to the K test. Compare against `besseli0`/`besseli1`, netlib `ribesl`, and tabulated reference.
- `test_besseli_negative_nu_reflection` — `besseli(-2.5_BK, x)` must agree with the closed-form reflection involving `besselk(2.5_BK, x)` to `eps(BK) * 10`.
- `test_besseli_integer_negative_x` — `besseli(3, -2.0_BK) == -besseli(3, 2.0_BK)` (odd `n`), `besseli(2, -2.0_BK) == besseli(2, 2.0_BK)` (even `n`).
- `test_besseli_non_integer_negative_x_nan` — `besseli(2.5_BK, -1.0_BK)` returns NaN.
- `test_besselix_no_overflow` — `besselix(0.0_BK, 700.0_BK)` finite while `besseli(0.0_BK, 700.0_BK)` overflows to `+Inf`.
- `test_besseli_nu_cputime`.

## Done when

- `besseli`, `besselix` exported.
- All four argument-domain branches (positive/negative `x`, integer/non-integer `nu`) covered with tests.
- Negative-`nu` reflection through `besselk` produces machine-precision agreement at integer-`nu` limits.
- Netlib `ribesl` benchmark shows >10x speedup at `nu=2.5, x=10` (current `besseli0` already beats netlib `ribesl` by ~200x; variable order will be slower per call but should still be dominant).
