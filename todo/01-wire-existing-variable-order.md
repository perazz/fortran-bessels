# 01 — Wire up existing variable-order machinery

## Goal

Ship `bessely(nu, x)`, `besselh(nu, k, x)`, `hankelh1(nu, x)`, `hankelh2(nu, x)`, and re-export `gamma_BK` from the `bessels` module. **No new numerical code** — only sign/parity wrappers around helpers that already exist.

## What's already there

Internal, private helpers in [src/bessels.f90](../src/bessels.f90):

| Helper | Line | What it does |
|---|---|---|
| `besselj_positive_args(nu, x)` | [252](../src/bessels.f90#L252) | Variable-order J for `nu ≥ 0`, `x ≥ 0`. All five branches: Debye, large-arg, Hankel-Debye, series, recurrence. |
| `bessely_positive_args(nu, x)` | [392](../src/bessels.f90#L392) | Variable-order Y for `nu ≥ 0`, `x > 0`. Integer-`nu` fast path + Debye + large-arg + Hankel + series + Chebyshev fallback. |
| `hankel_debye(nu, x) -> complex(BK)` | [src/bessels_debye.f90:57](../src/bessels_debye.f90#L57) | Returns `H^(1)_ν(x) = J_ν + i·Y_ν` via Debye expansion, valid for `x > ~nu`. |
| `bessely_power_series` | [src/bessels_constants.f90:762](../src/bessels_constants.f90#L762) | Returns `[Y_ν, J_ν]` simultaneously — useful for Hankel construction in the small-x branch. |

`besseljn(nu, x)` in [src/bessels.f90:183](../src/bessels.f90#L183) is the integer-order wrapper around `besselj_positive_args` and is already public — model `bessely(nu, x)` directly on its shape.

## What's missing

1. Public `bessely(nu, x)` wrapper (real `nu`, handles `x<0` → NaN, `x=0` → −∞, `nu<0` via reflection).
2. Public `besselh(nu, k, x)` wrapper returning `complex(BK)`, dispatching on `k ∈ {1, 2}`.
3. Public `hankelh1(nu, x)`, `hankelh2(nu, x)` — thin one-liners over `besselh`.
4. `bessels` module should `use bessels_gamma, only: gamma_BK` and re-export it so users get all advertised functions through one `use bessels`.

## Approach

### `bessely(nu, x)` — port [Bessels.jl/src/BesselFunctions/bessely.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/bessely.jl)

```fortran
elemental real(BK) function bessely(nu, x)
   real(BK), intent(in) :: nu, x
   real(BK) :: anu, Y, J
   integer  :: int_nu
   if (x < ZERO) then
       bessely = ieee_value(bessely, ieee_quiet_nan); return    ! Y undefined for x<0
   elseif (x == ZERO) then
       bessely = ieee_value(bessely, ieee_negative_inf); return
   end if
   anu = abs(nu)
   if (nu >= ZERO) then
       bessely = bessely_positive_args(nu, x)
   else
       ! Y_{-ν}(x) = cos(πν) Y_ν(x) + sin(πν) J_ν(x)
       Y = bessely_positive_args(anu, x)
       J = besselj_positive_args(anu, x)
       bessely = cos(PI*anu)*Y + sin(PI*anu)*J
   end if
end function
```

For integer `nu < 0`, `sin(πν) = 0` exactly (or near-zero), so the formula collapses to `Y_{-n} = (-1)^n Y_n` — the existing `cos(πν)` factor handles parity correctly if we let the elemental call short-circuit `nu ≥ 0` first. No special-case needed unless benchmarks show cancellation issues near integer `nu`.

### `besselh(nu, k, x)` — port [Bessels.jl/src/BesselFunctions/hankel.jl](https://github.com/heltonmc/Bessels.jl/blob/master/src/BesselFunctions/hankel.jl)

Two implementation paths:

- **`x > ~nu`** (Hankel-Debye regime): call `hankel_debye(nu, x)` directly — it already returns `J + i·Y`.
- **Otherwise**: compute `J = besselj_positive_args(nu, x)` and `Y = bessely_positive_args(nu, x)` separately, combine. (Bessels.jl has a fused `besseljy(nu,x)` — small optimization, defer.)

```fortran
elemental complex(BK) function besselh(nu, k, x)
   real(BK), intent(in) :: nu, x
   integer,  intent(in) :: k
   real(BK)    :: J, Y
   complex(BK) :: H
   if (hankel_debye_cutoff64(nu, x)) then
       H = hankel_debye(nu, x)
   else
       J = besselj_positive_args(nu, x)
       Y = bessely_positive_args(nu, x)
       H = cmplx(J, Y, BK)
   end if
   if (k == 1) then
       besselh = H
   else                 ! k == 2: H^(2) = conjg(H^(1)) for real x
       besselh = conjg(H)
   end if
end function

elemental complex(BK) function hankelh1(nu, x); ...; hankelh1 = besselh(nu, 1, x); end
elemental complex(BK) function hankelh2(nu, x); ...; hankelh2 = besselh(nu, 2, x); end
```

Negative-`nu` reflection for Hankels follows from `H^(k)_{-ν}(x) = exp(±iπν) H^(k)_ν(x)` — apply at the top of `besselh` before dispatching.

## Step-by-step

1. **Re-export `gamma_BK`.** In [src/bessels.f90](../src/bessels.f90) add `use bessels_gamma, only: gamma_BK` and `public :: gamma_BK`. Verify the existing test (`test_gamma` already calls it via `use bessels_gamma`) still passes; then update the test to `use bessels` only and confirm.
2. **Add `bessely(nu, x)`** as shown above, exported from `bessels`. Place near `besseljn` in [src/bessels.f90](../src/bessels.f90).
3. **Add `besselh`, `hankelh1`, `hankelh2`** after `bessely`. Export all three.
4. **Bench-guard**: confirm `bessely(0.0_BK, x)` and `bessely(1.0_BK, x)` match `bessely0(x)` / `bessely1(x)` to machine precision (they should — `bessely_positive_args` already special-cases integer `nu < 250` via `besselj_up_recurrence` from `bessely0/bessely1`).

## Tests (add to [test/bessels_test.f90](../test/bessels_test.f90))

- `test_bessely_nu` — for `nu ∈ {0.0, 0.5, 1.0, 1.5, 2.0, 3.7, 10.0}`, `x ∈ {0.1, 1.0, 5.0, 20.0, 100.0}`, compare against `bessely0` / `bessely1` for the integer-0/1 cases and against intrinsic `bessel_yn(int(nu), x)` for integer `nu`. Use Bessels.jl values from a known-good reference for non-integer `nu`.
- `test_hankelh1` — verify `hankelh1(0, x) = besselj0(x) + i*bessely0(x)` for `x ∈ {0.5, 5.0, 30.0}` covers all three Hankel branches.
- `test_hankelh2` — verify `hankelh2(nu, x) = conjg(hankelh1(nu, x))` for several `(nu, x)`.
- `test_bessely_nu_cputime` and `test_hankelh1_cputime` — add to the benchmark block at the bottom of the test file.

## Done when

- All four functions (`bessely`, `besselh`, `hankelh1`, `hankelh2`) plus `gamma_BK` are reachable via `use bessels` only.
- Test suite passes; new tests added.
- README "not yet implemented" list shrinks by four entries; "Currently available" list grows by four.
- No regression in existing ns/eval timings (these are new functions, not changes to hot paths).
