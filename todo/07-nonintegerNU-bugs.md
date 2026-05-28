# 07 — Pre-existing bugs in non-integer-ν paths

## Goal

Fix three port bugs in [src/bessels.f90](../src/bessels.f90), [src/bessels_constants.f90](../src/bessels_constants.f90), and [src/bessels_debye.f90](../src/bessels_debye.f90) that make `bessely(non-integer nu, x)` return wildly wrong values. The public `bessely(nu, x)` wrapper shipped in [01-wire-existing-variable-order.md](01-wire-existing-variable-order.md) is correct for integer ν only because of these bugs.

## Discovered when

Implementing plan 01 — added a half-integer reference test for `bessely(nu, x)` against the closed-form `Y_{1/2}/Y_{3/2}/Y_{5/2}` identities. Every single non-integer test case failed by 1-3 orders of magnitude. Walking through the call graph revealed three distinct port bugs.

The integer-ν path was unaffected (and still is) because it routes exclusively through `besselj_up_recurrence` from `bessely0/bessely1`, which doesn't touch any of the buggy helpers.

## Bug 1 — `bessely_power_series` argument order swap

[src/bessels_constants.f90:762](../src/bessels_constants.f90#L762):

```fortran
pure function bessely_power_series(x, nu) result(YJ)
   real(BK), intent(in) :: x, nu
```

[src/bessels.f90:424](../src/bessels.f90#L424):

```fortran
YJ = bessely_power_series(nu, x)   ! ← caller passes (nu, x), function signature is (x, nu)
```

The function body uses its local `x` as the argument (in `xo2 = HALF*x`) and local `nu` as the order (in `xo2**nu`). With the swap, the function computes `Y_{x_caller}(nu_caller)` instead of `Y_{nu_caller}(x_caller)`.

The Julia source has signature `bessely_power_series(v, x)` (v = nu first). The Fortran signature should be `(nu, x)` — but the call site is correct and the signature is wrong. Easier fix: swap the signature.

**Fix:** rename the parameters and/or swap in either the caller or the signature. Recommended:

```fortran
pure function bessely_power_series(nu, x) result(YJ)
   real(BK), intent(in) :: nu, x
```

(matches the Julia and the call site).

## Bug 2 — `bessely_chebyshev` mapping mismatch

[src/bessels.f90:415](../src/bessels.f90#L415):

```fortran
elemental complex(BK) function bessely_chebyshev(nu, x) result(cheb)
   ...
   nu_floor = nu - int(nu)
   Y = bessely_chebyshev_low_orders(nu_floor, x)
   call besselj_up_recurrence(x, Y%im, Y%re, nu_floor + ONE, nu, cheb%re, cheb%im)
```

`bessely_chebyshev_low_orders` returns a pair packed as `cheb%re, cheb%im` from `clenshaw(nu - 1, ...)` and `clenshaw(nu, ...)`. The recurrence then expects `(jnu at order nu_start-1, jnup1 at order nu_start)`. The current code passes `Y%im` first and `Y%re` second — these are at orders that don't match `nu_start = nu_floor + 1`.

This needs careful comparison to the Julia source (`src/BesselFunctions/bessely.jl` in Bessels.jl) to determine the correct argument ordering and the actual chebyshev table convention. Specifically:
- What range of `nu` does the table represent?
- Does `clenshaw_chebyshev(p, a)` evaluate at orderindex `p`, at `p + offset`, or at a normalized chebyshev point?
- In what order should the two table evaluations feed the recurrence?

Once that's pinned down, the fix is likely 1-2 lines. The recurrence itself is now correct (post-plan-06 bug fix) — only the call site needs adjustment.

## Bug 3 — `hankel_debye` complex-output truncation

[src/bessels_debye.f90:57-81](../src/bessels_debye.f90#L57):

```fortran
ab    = Uk_poly_Hankel(p, nu, -p2, x) ! why p*im ?
Uk_Yn = IM*ab(2)
hankel_debye = coef_Yn * Uk_Yn
```

The Julia source:

```julia
_, Uk_Yn = Uk_poly_Hankel(p*im, v, -p², T(x))
return coef_Yn * Uk_Yn
```

Two differences:
1. Julia passes `p*im` (complex first argument) to `Uk_poly_Hankel`. Fortran passes real `p`. The downstream `split_evalpoly` therefore can't reconstruct the correct complex polynomial value.
2. Even discarding the `*im` issue, the Fortran combines as `Uk_Yn = IM*ab(2)` (purely imaginary). Julia uses the second return directly, with `ab(2)` already a complex value from the complex polynomial evaluation.

For ν=0 inputs this also divides by zero (`-p/v = -0/0 = NaN`) in `Uk_poly{10,20}_split`, producing NaN outputs that leak through the test suite as `sum(z)=NaN`.

**Fix scope:** larger. `Uk_poly_Hankel` and the underlying `split_evalpoly` need a complex-argument variant. Either:
- Add a parallel `Uk_poly_Hankel_complex` returning `complex(BK)` for the Hankel-specific path, or
- Generalize the existing real path to accept and return complex via two real components, matching the Julia trick where (a, b) = real and imaginary parts of `poly(i*p)` evaluated by the even/odd split.

The Bessels.jl `split_evalpoly` in `Math.jl` shows the exact arithmetic — port that carefully.

**Workaround in shipped plan 01:** `besselh` was modified to bypass `hankel_debye` entirely (always uses `J + i·Y` composition). This makes `besselh(integer ν, x)` correct at all `x`, at the cost of speed in the large-x regime where `hankel_debye` was designed to be the fast path. `bessely_positive_args` still routes to `aimag(hankel_debye(...))` for non-integer ν in some x ranges — so non-integer `bessely(nu, x)` is still broken via this path too.

## Step-by-step

1. **Bug 1 first** — it's a one-liner swap. Verify by enabling the `test_bessely_nu_half_integer` test (removed in plan 01's restricted form) and confirming all `Y_{1/2}/Y_{3/2}/Y_{5/2}` cases pass at `x ∈ {0.5, 2.0}` (the power-series branch).
2. **Bug 2** — needs Julia-source comparison. Cross-check the table convention by computing a single chebyshev value by hand against `bessely(0.5, 7.0)` (a known value via closed form). Then fix the call.
3. **Bug 3** — port `Uk_poly_Hankel` with proper complex polynomial evaluation. Re-enable the `hankel_debye` fast path in `besselh` once verified. Add a `test_hankel_debye_correctness` test that compares `hankel_debye(nu, x)` to `cmplx(besselj_positive_args(nu, x), bessely_positive_args(nu, x))` for non-integer ν across the cutoff regions.
4. **Re-add the suppressed tests** in [test/bessels_test.f90](../test/bessels_test.f90): the half-integer ν test and a negative-ν Hankel reflection test (both removed in plan 01's restricted form to keep CI green).
5. **Restore the hankel_debye fast path** in `besselh` — search for the `! Implementation note:` comment block in [src/bessels.f90](../src/bessels.f90) and restore the `if (hankel_debye_cutoff(...)) then H = hankel_debye(...)` branch.

## Tests to re-enable / add

Removed from plan 01 to ship a clean suite:
- `test_bessely_nu_half_integer` — closed-form `Y_{1/2}/Y_{3/2}/Y_{5/2}`. Cover `x ∈ {0.5, 2.0, 7.0, 30.0, 100.0}` so each branch is exercised.
- `test_hankelh_negative_nu` — verify `H^(1)_{-ν}(x) = exp(+iπν)·H^(1)_ν(x)` for `ν ∈ {0.5, 2.5}`, `x = 5.0`.

New:
- `test_bessely_nu_chebyshev` — focused on `x ∈ (6, 19)` where the chebyshev branch is active.
- `test_hankel_debye_consistency` — direct check of `hankel_debye(nu, x)` vs `J + i·Y` for `(ν, x)` covering both Uk_poly10 and Uk_poly20 regimes.

## Done when

- All three bugs are fixed, with regression tests demonstrating each.
- The `besselh` workaround (bypass `hankel_debye`) is removed; the fast path is restored.
- `bessely(non-integer nu, x)` agrees with closed-form / tabulated reference to `eps(BK) * 100` across all branches.
- The suppressed tests are restored and pass.

## References

- Plan 01: [01-wire-existing-variable-order.md](01-wire-existing-variable-order.md) — where these bugs surfaced
- Bessels.jl reference: `src/BesselFunctions/bessely.jl`, `src/BesselFunctions/U_polynomials.jl`, `src/Math/Math.jl`
