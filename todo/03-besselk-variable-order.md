# 03 — `besselk(nu, x)` variable order

## Goal

Public `besselk(nu, x)` and scaled `besselkx(nu, x) = exp(x) · K_ν(x)`, both `elemental real(BK)`, real `nu` and `x`.

This is the **largest** of the variable-order ports: it introduces two pieces of new numerical infrastructure (Temme series for K near integer order, and the K-specific Debye uniform asymptotic) that don't have analogues in the existing code.

## What's already there

- `besselk0(x)`, `besselk1(x)` in [src/bessels.f90](../src/bessels.f90) for the base cases of forward recurrence.
- `gamma_BK` in [src/bessels_gamma.f90](../src/bessels_gamma.f90) for the Temme series.
- `evalpoly`, `muladd`, `cbrt` utilities.
- The `Uk_poly{5,10,20}` machinery in [src/bessels_debye.f90](../src/bessels_debye.f90) — note that these are the J/Y Debye polynomials. K_ν uses a **different** set of U-coefficients (NIST 10.41); we need to port those too. They share the same recursive structure but the leading-order behavior differs.

## What's missing

1. K-specific Debye uniform asymptotic expansion — call it `besselk_debye(nu, x)`. Distinct from `besseljy_debye`; see NIST DLMF 10.41 and Bessels.jl `src/BesselFunctions/U_polynomials.jl` (`Uk_poly_Kn`).
2. Temme series for K near integer order — Bessels.jl `src/BesselFunctions/besselk.jl` `_K_nu_temme`. Handles the `0/0` cancellation when `ν → n` exactly via L'Hôpital-style limit.
3. Large-x asymptotic `K_ν(x) ≈ √(π/(2x)) · e^{-x} · P(1/x)` (asymptotic series in `1/x`, coefficients functions of `ν`).
4. Forward recurrence wrapper — `K_{n+1} = (2n/x) K_n + K_{n-1}` (note `+`, opposite sign from J's recurrence — both directions are stable for K because K grows downward).
5. Public `besselk(nu, x)` dispatching across the four branches.
6. Scaled `besselkx(nu, x)`.

## Approach — port [Bessels.jl/src/BesselFunctions/besselk.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/besselk.jl)

Four branches by `(nu, x)`:

| Branch | Method |
|---|---|
| `x > nu²/36 + 18` | Large-x asymptotic series in `1/x`. |
| `x > nu⁴/2401 + 1.5` (and not branch 1) | Levin-transform sequence acceleration on the asymptotic series. **Defer this — the other three branches cover most of parameter space; mark as a follow-up.** |
| `nu > 25` or `x > 35` (and not above) | Debye uniform asymptotic via K-specific U-polynomials. |
| Default (small `x`, moderate `nu`) | Temme series for non-integer `nu`; for integer `nu` use upward recurrence from `besselk0`/`besselk1`. |

`K_{-ν} = K_ν` exactly — handle at the top with `nu = abs(nu)`, no reflection cost.

`x ≤ 0`: NaN (matches `besselk0`'s convention).

## Step-by-step

1. **Negative-`nu` and edge cases.** Add the top-level dispatcher returning NaN for `x ≤ 0` and folding `nu = abs(nu)`.
2. **Integer-`nu` fast path.** Use the existing `besselk0`/`besselk1` and upward recurrence `K_{n+1}(x) = (2n/x) K_n(x) + K_{n-1}(x)`. Add a helper `besselk_up_recurrence` mirroring `besselj_up_recurrence` in [src/bessels_constants.f90:843](../src/bessels_constants.f90#L843) — same shape, plus instead of minus in the recurrence. **Watch out** for the suspected bug in `besselj_up_recurrence` (see [06-housekeeping.md](06-housekeeping.md)) — don't propagate it.
3. **Large-x asymptotic.** Port `K_large_argument` from Bessels.jl. The leading factor is `√(π/(2x)) · e^{-x}`; the series in `1/x` has coefficients `μ = 4ν²`, then `a_k = ∏(μ - (2j-1)²)/(k! · 8^k · x^k)`. For `besselkx`, drop the `e^{-x}` factor. Use Horner via `evalpoly`; the cutoff is when adding the next term would not change the sum to `eps(BK)`.
4. **Temme series.** This is the non-trivial new component. Port from `_K_nu_temme` — it computes `K_ν` and `K_{ν+1}` for non-integer `ν` near zero, then can be upward-recurred. The series converges for `x ≤ 2`, and combines with downward recurrence for larger `x`. Key trick: a power-series representation of `K_ν` that stays accurate even as `ν → integer` (the singularity in `gamma(-ν)` cancels analytically).
5. **K-specific Debye expansion.** Port `besselk_debye` (and its U-polynomial table) from `U_polynomials.jl`. The same `Uk_poly_{5,10,20}` shape applies but with the K-flavored coefficients. **Plan to share infrastructure with [04-besseli-variable-order.md](04-besseli-variable-order.md)** — `besseli_debye` uses the *same* U-table (just a different prefactor). Put the new K-flavored U-polys in [src/bessels_debye.f90](../src/bessels_debye.f90) next to the existing J/Y ones, or split into `bessels_debye_ki.f90` if the file gets unwieldy.
6. **Dispatcher.** Branch as per the table above. Add cutoff helpers `besselk_debye_cutoff`, `besselk_large_arg_cutoff` mirroring the existing `besseljy_debye_cutoff64` etc.
7. **`besselkx`.** Same dispatcher; in each branch skip or add back the `e^{±x}` factor as appropriate.

## New files / file changes

- [src/bessels_besselk.f90](../src/bessels_besselk.f90) — new module hosting `besselk_debye`, `besselk_temme_series`, `besselk_large_argument`, and the public `besselk` / `besselkx`. Keeps [src/bessels.f90](../src/bessels.f90) from ballooning.
- [src/bessels_debye.f90](../src/bessels_debye.f90) — add K U-polynomial tables.
- [src/bessels.f90](../src/bessels.f90) — `use bessels_besselk`, export `besselk`, `besselkx`.

## Tests

- `test_besselk_nu` — values at `(nu, x)` grid `nu ∈ {0.0, 0.5, 1.0, 2.5, 5.0, 30.0, 100.0}`, `x ∈ {0.01, 0.1, 1.0, 5.0, 20.0, 100.0, 500.0}`. Compare against `besselk0`/`besselk1` for `nu=0/1`, against netlib `rkbesl` (already in [test/3rd_party/](../test/3rd_party/)) for integer `nu`, and against tabulated reference for non-integer.
- `test_besselk_temme` — focused on `x < 2`, `nu ∈ {0.1, 0.4, 0.99, 1.01, 2.0001}` — verify the near-integer limit doesn't blow up.
- `test_besselk_debye` — `nu = 100`, `x ∈ {5, 50}` — verify the Debye branch is active and accurate.
- `test_besselkx_overflow` — `besselkx(0.0, 700.0)` should be finite (~`√(π/1400)`) while `besselk(0.0, 700.0)` underflows to zero.
- `test_besselk_nu_cputime` benchmark.

## Done when

- `besselk(nu, x)` and `besselkx(nu, x)` exported from `bessels`.
- Tests pass; netlib comparison shows we're at least as fast (current `besselk0/1` are ~5-6 ns/eval — variable-order will be slower, target `<30 ns/eval` for `besselk_debye` branch).
- The K U-polynomial machinery is shared with the upcoming `besseli` port.
- Levin acceleration deferred to a follow-up TODO; document the cutoff gap (region between "large-x asymptotic" and "Debye" where accuracy may be ~1 ulp worse).
