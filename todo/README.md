# Implementation Roadmap

Prioritized plans for the functions still missing from the public API, based on what's already implemented internally and what dependencies each new function brings.

## State of the codebase (audit, 2026-05-27)

**Public API today** — [src/bessels.f90](../src/bessels.f90) exports `besselj0/j1/jn`, `bessely0/y1`, `besselk0/k1`, `besseli0/i1`, plus `cbrt` / constants. `gamma_BK` lives in [src/bessels_gamma.f90](../src/bessels_gamma.f90) but is **not** re-exported from `bessels` — users have to `use bessels_gamma` separately, which contradicts the README.

**Hidden assets** — already implemented as private helpers but not exported:
- `besselj_positive_args(nu, x)`, `bessely_positive_args(nu, x)` — full variable-order J/Y machinery (debye + large-arg + hankel + series + recurrence branches).
- `hankel_debye(nu, x) -> complex(BK)` — full Hankel implementation for `x > ~nu`.
- All `Uk_poly{5,10,20}` Debye U-polynomial tables, `besseljy_debye`, `besseljy_large_argument`, `clenshaw_chebyshev`, `bessely_chebyshev*`, power series for J/Y, up/down recurrences.

This means **Tier 0** below is mostly plumbing.

## Suggested order

| #  | Plan | Effort | Why this order |
|----|------|--------|----------------|
| 0  | [01-wire-existing-variable-order.md](01-wire-existing-variable-order.md) | small | Ship `bessely(nu,x)`, `besselh`, `hankelh1/2`, re-export `gamma_BK` — most of the work is already done. |
| 1  | [02-airy.md](02-airy.md) | medium | Fully independent of Bessel routines (Bessels.jl uses Cephes-style minimax polys, not K_{1/3}). Can land anytime. |
| 2  | [03-besselk-variable-order.md](03-besselk-variable-order.md) | large | New infrastructure: Temme series + uniform asymptotic for K_ν. `besseli(nu,x)` depends on it for negative-ν reflection. |
| 3  | [04-besseli-variable-order.md](04-besseli-variable-order.md) | medium | Depends on K_ν. Power series + Debye + large-x asymptotic + scaled variant `besselix`. |
| 4  | [05-spherical-bessels.md](05-spherical-bessels.md) | medium | Needs Tier 0–3 done. Mix of integer-order fast paths (sin/cos/sinh/cosh closed forms) and half-integer fallback. |
| 5  | [06-housekeeping.md](06-housekeeping.md) | small | Four in-code TODOs in `bessely0`/`bessely1`, suspected bug in `besselj_up_recurrence`, `Uk_poly_Hankel` "why p*im?" comment. |
| 6  | [07-nonintegerNU-bugs.md](07-nonintegerNU-bugs.md) | medium | Three port bugs (`bessely_power_series` arg swap, `bessely_chebyshev` mapping, `hankel_debye` complex output) surfaced while implementing plan 01. `bessely(nu, x)` and `besselh` ship correct for integer ν only until this lands. |

## Out of scope (matching Bessels.jl's public API)

- Scaled variants `besseli0x`, `besselix`, `besselk0x`, `besselkx`, `airyaix`, etc. — fold into each tier as the unscaled version lands.
- In-place sequence variants (`besselj!`, `bessely!`, …) — useful for HPC callers; revisit after the scalar API is complete.
- Complex Airy (`cairy.jl`) and complex-argument ordinary Bessels — **not** in Bessels.jl's stable public API; defer indefinitely.

## How each plan is structured

Each file follows the same shape:

1. **Goal** — one sentence, what the public-facing change is.
2. **What's already there** — concrete file/function references for code we can reuse.
3. **What's missing** — what we have to write or port.
4. **Approach** — branch structure, citing the Julia reference file.
5. **Step-by-step** — numbered, each step verifiable.
6. **Tests** — what to add to `test/bessels_test.f90`.
7. **Done when** — measurable success criteria.

Reference: Bessels.jl at https://github.com/heltonmc/Bessels.jl (`src/BesselFunctions/`, `src/AiryFunctions/`, `src/Math/`).
