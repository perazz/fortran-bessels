# 06 — Housekeeping & suspected bugs

Small fixes, leftover TODOs, and one bug found during the audit. None of these block the variable-order work, but they should land alongside it.

## 1. `besselj_up_recurrence` — suspected bug

[src/bessels_constants.f90:843-871](../src/bessels_constants.f90#L843)

```fortran
elemental subroutine besselj_up_recurrence(x, jnu, jnum1, nu_start, nu_end, a, b)
    ...
    nu   = nu_start
    do while (nu<nu_end+HALF)
        jnu2 = [jnu2(2), nu_start*x2*jnu2(2) - jnu2(1)]   ! ← nu_start, not nu
        nu = nu-ONE                                          ! ← decrement, not increment
    end do
```

Forward recurrence for J/Y is `f_{n+1} = (2n/x)·f_n − f_{n−1}`. The coefficient should be `2n/x = nu * x2` where `nu` increases each iteration — but the code uses the *fixed* `nu_start` and *decrements* `nu`. That makes the `while (nu<nu_end+HALF)` condition only true if `nu_start < nu_end`, in which case the loop never advances `nu` toward termination (it goes the wrong way) → infinite loop, **or** the loop runs once and returns the wrong value.

**Why it's not blowing up today**: `besselj_up_recurrence` is only called from `bessely_positive_args` (integer `nu` fast path) and `bessely_chebyshev`. The integer-`nu` path calls with `nu_start = ONE, nu_end = nu` — when `nu == 1` the loop body runs once (with the wrong coefficient, but for a degenerate case it produces the right answer) and exits. For `nu > 1` the test suite would catch it — except `test_bessel_y0`/`y1` don't actually exercise variable-order Y, and `bessely(nu, x)` isn't public yet, so the bug is dormant.

**Fix:** replace `nu_start` with `nu` inside the loop and increment `nu = nu + ONE`. Add a regression test that compares `besselj_up_recurrence(x, j1, j0, 1.0, 5.0, ...)` against five forward applications of the recurrence by hand.

**Cross-check the down recurrence**: [src/bessels_constants.f90:822-837](../src/bessels_constants.f90#L822) `besselj_down_recurrence` uses `nus*x2*jnu2(2)-jnu2(1)` with `nus = nu_start` then `nus = nus-ONE` — but it's down-recurrence (`nu_start > nu_end`) so the coefficient should decrease each step, which `nus = nus - ONE` does. The down version looks correct; only the up version is broken.

This must be fixed **before** [01-wire-existing-variable-order.md](01-wire-existing-variable-order.md) ships, because `bessely(nu, x)` will route integer-`nu` calls through this broken recurrence.

## 2. Four in-code `TODO` comments

[src/bessels.f90](../src/bessels.f90):

- Line 488: `! TODO: replace the two polynomials with a single one` (in `bessely0`, small-x branch)
- Line 497: same TODO in `bessely0`, mid-x branch
- Line 560: same TODO in `bessely1`
- Line 570: same TODO in `bessely1`

These are about rational-function `P/Q` evaluations where `P` and `Q` are evaluated independently then divided. The TODO suggests fusing into a single evaluation. **Defer** unless benchmarks show `bessely0/y1` are hot in caller workloads — they're already at 17-18 ns/eval (faster than the gfortran intrinsic). Convert to GitHub issues or strike out if not pursuing.

## 3. `gamma_BK` not re-exported

[src/bessels.f90](../src/bessels.f90) does not `use bessels_gamma`. The README implies `gamma_BK` is part of the package's API, and the test file does `use bessels_gamma` directly to get at it. This split is awkward. Pick one:

- **Re-export from `bessels`** (recommended): add `use bessels_gamma, only: gamma_BK` and `public :: gamma_BK`. One-line change, no behavioral risk. Bundled with [01-wire-existing-variable-order.md](01-wire-existing-variable-order.md).
- **Or split clearly** — document that math helpers live in `bessels_gamma` and update the README accordingly.

Going with the re-export keeps `use bessels` as the single entry point.

## 4. `Uk_poly_Hankel` — leftover comment

[src/bessels_debye.f90:76](../src/bessels_debye.f90#L76):

```fortran
ab    = Uk_poly_Hankel(p, nu, -p2, x) ! why p*im ?
```

Someone (Federico?) left a confused comment while porting. `hankel_debye` produces a complex output and the U-polynomial values are split into the real/imaginary parts via `split_evalpoly` — the `-p2` argument is correct (the alternating-sign split). The `p*im` reference in the comment is a Julia-ism that didn't translate; the Fortran path computes `Uk_Yn = IM * ab(2)` on line 77 which is the correct construction. **Delete the comment.**

## 5. Multi-precision cutoffs declared but unused

[src/bessels_constants.f90](../src/bessels_constants.f90) has `besseljy_debye_cutoff32`, `besseljy_debye_cutoff128`, `bessely_series_cutoff32`, `besselj_series_cutoff128` etc. — `real32` and `real128` variants of every cutoff. None are wired up because `BK` is fixed to `real64`. Either:

- Keep them (low cost, signals future multi-precision intent — the `! Todo: make one module per real precision` comment at [src/bessels.f90:24](../src/bessels.f90#L24) confirms the intent).
- Move to a `bessels_constants_multi.f90` to keep the `real64` module lean.

**Defer**. Revisit when a multi-precision branch actually starts.

## Step-by-step

1. Fix `besselj_up_recurrence` (item 1). Add the regression test. Run the existing test suite to confirm nothing else regresses.
2. Add `use bessels_gamma, only: gamma_BK` and `public :: gamma_BK` to [src/bessels.f90](../src/bessels.f90) (item 3). Update the test file to `use bessels` only for `gamma_BK`.
3. Delete the stray `! why p*im ?` comment (item 4).
4. Decide on items 2 and 5; either fix or close as won't-fix and remove the inline TODOs.

## Done when

- `besselj_up_recurrence` produces correct values for `nu > 1`, verified by direct comparison.
- `gamma_BK` reachable via `use bessels`.
- Cleanup PR is small and surgical — no behavioral changes to anything other than `besselj_up_recurrence`.
