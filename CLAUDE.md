[//]: # (CLAUDE.md — fortran-bessels)

# CLAUDE.md — fortran-bessels

An open-source, MIT-licensed modern Fortran port of [Bessels.jl](https://github.com/heltonmc/Bessels.jl). The aim is to be **as accurate as possible and as fast as possible** in pure Fortran, using modern Fortran style. Performance is a first-class concern; every change should preserve or improve both accuracy and speed.

## Core Principles

**1. Speed & accuracy are the goal.** The library should beat the compiler intrinsics and netlib `specfun` on the benchmarks in [README.md](README.md). Don't trade measurable performance for marginal clarity; do trade marginal performance for measurable accuracy.
**2. Simplicity first.** Minimum code that solves the problem. No speculative abstractions, no configurability for hypothetical callers, no error handling for impossible cases. Branch-on-`x` polynomial evaluation is the dominant pattern — don't dress it up.
**3. Surgical changes.** Touch only what the task requires. Match the existing style (banner header, `elemental real(BK) function`, branch comments before the function body). Don't reformat adjacent code.
**4. Match Bessels.jl semantics.** When in doubt about a branch cutoff, an asymptotic expansion, or a coefficient table, the Julia source is the reference. Cite the Julia file or paper in a one-line comment when reproducing a non-obvious formula.
**5. Verify with the test suite.** `bessels_test.f90` covers correctness (against intrinsics / netlib) and prints ns/eval timings. Any change that touches a public function must run cleanly and not regress timings.

## Project Rules

- **Public GitHub repo** at `perazz/fortran-bessels`, MIT license. Contributions welcome via PR.
- **Commits and PRs** on the owner's account — no `Co-Authored-By` lines, no "Generated with Claude Code" attributions.
- **No version bumps per PR.** Release tags only.
- **Pure Fortran only.** The `test/3rd_party/` `.f90` files (netlib `ribesl`, `rkbesl`) are *test-time reference implementations* only — never call them from `src/bessels*.f90`.

## Building & Testing

The canonical build is [fpm](https://fpm.fortran-lang.org):

```bash
fpm test --profile release                 # build + run bessels_tests
fpm test --profile release --flag "-march=native"
```

Manual gfortran build (matches the README recipe):

```bash
gfortran -ffree-line-length-none -O3 -march=native -ffast-math \
    src/bessels_constants.f90 src/bessels_gamma.f90 src/bessels_debye.f90 \
    src/bessels.f90 \
    test/3rd_party/ribesl.f90 test/3rd_party/rkbesl.f90 \
    test/bessels_test.f90 -o bessels_test
./bessels_test
```

A Code::Blocks workspace exists at [project/fortran-bessels.cbp](project/fortran-bessels.cbp) and is kept in sync with the fpm layout.

**Performance flags matter.** Benchmark numbers in the README assume `-O3 -march=native -ffast-math`. When investigating a perf regression, always check with those flags — `-O0` or `-O2` will mislead.

**Test output policy.** `bessels_test.f90` prints a one-line pass/fail per test plus ns/eval for each benchmarked function. Tests should not print extra noise; gate any diagnostic prints behind a local `logical, parameter :: debug = .false.`.

## File & Module Layout

```
src/
  bessels_constants.f90   — kinds (BK, BSIZE), named constants, polynomial coefficient tables
  bessels_gamma.f90       — gamma helpers used by Bessel branches
  bessels_debye.f90       — Debye / uniform asymptotic expansions for large nu
  bessels.f90             — public API: besselj0/j1/jn, bessely0/y1, besseli0/i1, besselk0/k1
test/
  bessels_test.f90        — correctness + timing harness
  3rd_party/              — netlib reference (test only, do not call from src)
                            * ribesl.f90 → module bessels_ribesl
                            * rkbesl.f90 → module bessels_rkbesl
```

One module per file, each with the project's banner header:

```fortran
!  ************************************************************************************************************
!                                    ____  ________________ ________   _____
!                                   / __ )/ ____/ ___/ ___// ____/ /  / ___/
!                                  ...
!                                              <module purpose>
!  MIT License
!  Copyright (c) 2022-2026 Federico Perini
!  Copyright (c) 2021-2022 Michael Helton, Oscar Smith, and the Bessels.jl contributors
!  ************************************************************************************************************
```

## Style — Types, Constants, Procedures

**Kinds.** All real work is in `real(BK)`, where `BK = real64`, and all integer counters that need a specific width use `BSIZE = int32`. **Never** plain `real ::` / `double precision ::`. Define new constants with the `_BK` suffix: `1.0_BK`, `26.0_BK`.

**Named constants over literals.** The `bessels_constants` module exports `ZERO`, `ONE`, `HALF`, `THIRD`, `PI`, `PIO2`, `PIO4`, `TWOOPI`, `SQ2OPI`, `SQ1O2PI`, etc. Use them. Hardcoded magic literals are acceptable *only* when they're polynomial coefficients in a `parameter` array (where the value itself is the point).

**Procedure attributes.**
- `elemental real(BK) function f(x)` is the default shape for every public Bessel function. This is what lets callers pass arrays directly.
- `pure` for any helper without side effects.
- `intent(in)` / `intent(out)` / `intent(inout)` on every argument.

**Function skeleton** (matches `besselj0`, `besselj1`, `besseli0`, …):

```fortran
elemental real(BK) function besselXY(x)
   real(BK), intent(in) :: x
   real(BK) :: ax, xinv, x2, ...
   integer  :: n
   real(BK), parameter :: ppoly(*) = [ ONE, -1.0_BK/16.0_BK, ... ]
   ax = abs(x)
   if (ax <= CUTOFF1) then
       besselXY = evalpoly(size(POLY_SMALL), ax**2, POLY_SMALL)
   elseif (ax < CUTOFF2) then
       ! Branch comment: which expansion, around which root, why
       ...
   else
       ! Large-x asymptotic
       ...
   end if
end function
```

A short branch comment block immediately *before* the function body (see `besselj0` in [src/bessels.f90](src/bessels.f90)) is the established documentation pattern. Keep it.

**Use `iso_fortran_env`.** `use iso_fortran_env, only: real64, int32, ...` at the top of `bessels_constants`. Other modules pull kinds from `bessels_constants`.

## Control Flow, Arrays, Performance

**Branch on `abs(x)`.** Every Bessel function is structured as a piecewise polynomial / asymptotic expansion over disjoint ranges of `|x|`. Keep branches ordered from cheapest/most-common to expensive, and use `elseif` chains — not `select case` on reals (illegal) and not nested `if/then/else` trees.

**Polynomial evaluation.** Use the `evalpolyN` / `evalpoly` helpers and `muladd`. These are tuned for FMA contraction; do not hand-roll Horner schemes.

**Array intrinsics > explicit loops** where the function isn't already `elemental`. The elemental Bessel functions vectorize automatically over array arguments.

**Fortran `.and.` / `.or.` are non-short-circuiting** — never write `if (present(x) .and. x > 0)`. Always nest: `if (present(x)) then; if (x > 0) ...; end if`.

**Column-major layout.** For the coefficient tables in `bessels_constants` (`J0_POLYS(:,n)`, `J0_ROOTS(:,n)`), the *first* index is the inner loop. Don't transpose these tables.

## Numeric Edge Cases

Each public function must do the right thing for:
- `x = 0` (exact: `j0(0)=1`, `j1(0)=0`, `y0(0) = -inf`, `y1(0) = -inf`, `k0(0)=+inf`, `i0(0)=1`, …)
- Negative `x` for the J-family (`j0` even, `j1` odd, `jn` parity follows `n`)
- Negative `x` for Y/K/I where the function is undefined → return `ieee_quiet_nan`
- `x = +infinity` (asymptote to 0 for J/Y/K, to ±∞ for I)
- `x > huge(x)` and subnormals

Don't add bounds-checking branches that aren't load-bearing for accuracy — the test suite is the spec.

## GCC Pitfalls (gfortran 12+)

1. **Fortran `.and.`/`.or.` do not short-circuit.** Use nested `if`.
2. **`-ffast-math` reorders FMAs.** If a branch's accuracy depends on a specific evaluation order, isolate it in a `pure` helper and document the constraint in a one-line comment. Otherwise let the compiler reorder.
3. **`elemental` + branchy bodies vectorize poorly under some gfortran versions.** If a profiler shows a hot function not vectorizing, check `-fopt-info-vec-missed` before restructuring.
4. **Avoid `select case` on real-valued cutoffs** (illegal in standard Fortran).
5. **Watch parameter-array constructors.** `[ONE, -1.0_BK/16.0_BK, ...]` is fine for `real(BK), parameter`; do not introduce derived-type array constructors here.

## Do Nots

- Default kinds (`real ::`, `double precision ::`) — always `real(BK)`.
- Bare magic literals (`2.d0`, `3.14...`) outside polynomial coefficient tables — use `bessels_constants`.
- Call netlib `ribesl` / `rkbesl` from `src/` — they exist for `test/` reference only.
- Heavy abstractions (classes, type-bound dispatch) for the public API. The Bessel functions are free `elemental` functions — keep them that way.
- Performance regressions without an accompanying accuracy gain that justifies them. If a change costs ns/eval, the commit message must say why.
- `Co-Authored-By` / "Generated with Claude Code" lines in commits or PRs.
- Version bumps on every PR — release tags only.
