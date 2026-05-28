# 02 — Airy functions

## Goal

Public `airyai(x)`, `airyaiprime(x)`, `airybi(x)`, `airybiprime(x)` plus scaled variants `airyaix`, `airyaiprimex`, `airybix`, `airybiprimex` — all `elemental real(BK)`. Real argument only; complex Airy is out of scope.

## What's already there

Nothing. Airy lives in its own world: Bessels.jl deliberately **avoids** routing Airy through `K_{1/3}` / `J_{1/3}` (the textbook approach) for speed and accuracy reasons. So none of the existing Bessel infrastructure is reused.

Reusable from existing code:
- `evalpoly`, `evalpoly4/5/7/8`, `muladd`, `clenshaw_chebyshev` in [src/bessels_constants.f90](../src/bessels_constants.f90).
- `cbrt` for `x^{2/3}` operations: `x**(2.0_BK/3.0_BK)` or `cbrt(x)**2`.

## What's missing

Everything — a new module [src/bessels_airy.f90](../src/bessels_airy.f90), wired into [src/bessels.f90](../src/bessels.f90) via `use bessels_airy` and re-export.

## Approach — port [Bessels.jl/src/AiryFunctions/airy.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/AiryFunctions/airy.jl)

Four argument regions, the same for `Ai` and `Bi` (with different coefficient tables):

| Region | Method |
|---|---|
| `x ∈ [0, 2.06)` | Minimax polynomial / rational approximation (Cephes-derived). |
| `x ≥ 2.06` | Tabulated asymptotic polynomial expansion: `Ai(x) ≈ exp(-ζ)/(2√π x^{1/4}) · P(1/ζ)` with `ζ = (2/3)x^{3/2}`. Bessels.jl pre-computed `P` coefficients in Mathematica/Arb. |
| `x ∈ [-10, 0)` | Piecewise Taylor series around the Airy zeros. |
| `x < -10` | Asymptotic with Euler-formula trick to stay in real arithmetic: `Ai(x) ≈ (1/√π)·|x|^{-1/4}·[cos(ζ+π/4)·P + sin(ζ+π/4)·Q]` where `ζ = (2/3)|x|^{3/2}`. |

`Bi` uses the same branch structure but with sign flips and different coefficient tables; some intervals use `clenshaw_chebyshev`.

**Derivatives** (`airyaiprime`, `airybiprime`) have their own coefficient tables — they are not numerical derivatives of the value functions. Bessels.jl computes them with the same branch structure.

**Scaled variants** apply `exp(2/3 · x^{3/2})` (for Ai) or `exp(-2/3 · x^{3/2})` (for Bi). They **throw `DomainError` for negative x** in Bessels.jl — in Fortran return `ieee_quiet_nan` to match the existing convention (see `bessely0` for `x<0`).

## Step-by-step

1. **Create [src/bessels_airy.f90](../src/bessels_airy.f90).** Banner header, `use bessels_constants`, `private`, `public :: airyai, airyaiprime, airybi, airybiprime, airyaix, airyaiprimex, airybix, airybiprimex`.
2. **Port the coefficient tables.** From Bessels.jl `src/AiryFunctions/airy.jl` — copy as `real(BK), parameter` arrays. These are tabulated and stable; mass-port them. Group by function (`AI_SMALL_P`, `AI_SMALL_Q`, `AI_LARGE_P`, `AI_NEG_P`, etc.) following the Julia naming.
3. **`airyai(x)`** with the four branches. Each branch: a polynomial / rational call + a multiplicative prefactor. Use `cbrt` only where the cubic-root form is more accurate than `x**(1.0_BK/3.0_BK)` — they should be equivalent at `real64` but `cbrt` is what Bessels.jl uses for `x^{1/3}` integers-of-three exponents.
4. **`airybi(x)`** — same structure, different tables. Watch the `clenshaw_chebyshev` intervals: pass the Chebyshev coefficients as the second arg to `clenshaw_chebyshev`.
5. **`airyaiprime(x)`, `airybiprime(x)`** — separate coefficient sets, same branch structure.
6. **Scaled variants** — pass `do_scale` as a boolean parameter? Bessels.jl has them as separate top-level functions; mirror that. For each, branch on `x < 0` → NaN, else compute as `airyai(x) * exp(2x^{3/2}/3)` for Ai (and analogous for Bi). The scaled variants exist because the unscaled `Ai(x)` underflows for `x > ~104` and `Bi(x)` overflows; the scaled versions stay finite.
7. **Wire into [src/bessels.f90](../src/bessels.f90)**: `use bessels_airy`, add the eight functions to the `public` list.
8. **Edge cases.** `x = 0`: `Ai(0) = 1/(3^{2/3} Γ(2/3))`, `Bi(0) = 1/(3^{1/6} Γ(2/3))`, `Ai'(0) = -1/(3^{1/3} Γ(1/3))`, `Bi'(0) = 3^{1/6}/Γ(1/3)`. The constants module already has `GAMMA_TWO_THIRDS` and `GAMMA_ONE_THIRD`. Bessels.jl returns these exactly without branching — verify the small-`x` polynomial reproduces them to within `eps()`.

## Tests (add to [test/bessels_test.f90](../test/bessels_test.f90))

- `test_airyai` — values at `x ∈ {-50, -10, -2, -0.5, 0, 0.5, 2, 10, 50}` against tabulated reference (Wolfram Alpha / Bessels.jl).
- `test_airyaiprime`, `test_airybi`, `test_airybiprime` — same grid.
- `test_airyai_zeros` — `Ai(a_k) = 0` for the first few Airy zeros (`a_1 ≈ -2.3381…`, `a_2 ≈ -4.0879…`); verify the negative-x branch is accurate near them.
- `test_airyaix_overflow_guard` — `airyaix(100.0)` should be `O(1)` while `airyai(100.0)` underflows.
- `test_airyai_cputime`, `test_airybi_cputime` — benchmark block.

No netlib comparison — netlib `specfun` Airy is in `src/3rd_party/` only if we add it, and Bessels.jl avoids that route deliberately. Use tabulated reference values.

## Done when

- Eight Airy functions exported from `bessels`.
- All test cases pass to `eps(BK) * 10` relative error on the value functions, `eps(BK) * 100` near zeros (zeros are inherently ill-conditioned).
- Benchmarks added; ns/eval competitive with the rest of the library (Airy is polynomial-heavy → should be in the 20-40 ns/eval range).
- README "Not yet implemented" loses four entries; "Currently available" gains eight (four + four scaled).
