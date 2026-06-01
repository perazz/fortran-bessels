#!/usr/bin/env python3
"""Generate documentation plots for fortran-bessels.

Produces PNG images in doc/media/ for the Doxygen function pages:

  * Tier A — function-shape curves drawn with scipy.special (one per family).
  * Tier B — true accuracy plots: relative error (in ULPs of real64) of this
             library's own output against an mpmath arbitrary-precision
             reference. The library values are harvested by running the Fortran
             helper `example/dump_values.f90` via fpm.

Run from the repository root:

    python doc/generate_plots.py

Requires: numpy, matplotlib, scipy, mpmath  (Tier B also needs an fpm/gfortran
toolchain; if fpm is missing the accuracy plots are skipped).
"""
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

ROOT = Path(__file__).resolve().parent.parent
MEDIA = Path(__file__).parent / "media"
DATA = MEDIA / "data"
MEDIA.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)

EPS = np.finfo(np.float64).eps  # 1 ULP at 1.0

# Shared style — identical to the sibling fitpack doc site for a consistent look.
plt.rcParams.update({
    "figure.figsize": (6, 4),
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
})


# ---------------------------------------------------------------------------
# Tier A — function-shape curves (scipy.special)
# ---------------------------------------------------------------------------

def plot_besselj():
    x = np.linspace(0, 20, 800)
    fig, ax = plt.subplots()
    ax.plot(x, sp.j0(x), label="$J_0(x)$")
    ax.plot(x, sp.j1(x), label="$J_1(x)$")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlabel("x"); ax.set_ylabel("$J_\\nu(x)$")
    ax.set_title("Bessel functions of the first kind")
    ax.legend(); _save(fig, "besselj.png")


def plot_besseljn():
    x = np.linspace(0, 20, 800)
    fig, ax = plt.subplots()
    for n in range(5):
        ax.plot(x, sp.jv(n, x), label=f"$J_{n}(x)$")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlabel("x"); ax.set_ylabel("$J_n(x)$")
    ax.set_title("Integer-order Bessel functions $J_n$")
    ax.legend(ncol=2, fontsize=9); _save(fig, "besseljn.png")


def plot_bessely():
    x = np.linspace(0.05, 20, 800)
    fig, ax = plt.subplots()
    ax.plot(x, sp.y0(x), label="$Y_0(x)$")
    ax.plot(x, sp.y1(x), label="$Y_1(x)$")
    ax.set_ylim(-2.0, 0.8); ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlabel("x"); ax.set_ylabel("$Y_\\nu(x)$")
    ax.set_title("Bessel functions of the second kind")
    ax.legend(); _save(fig, "bessely.png")


def plot_besseli():
    x = np.linspace(0, 5, 600)
    fig, ax = plt.subplots()
    ax.plot(x, sp.i0(x), label="$I_0(x)$")
    ax.plot(x, sp.i1(x), label="$I_1(x)$")
    ax.set_xlabel("x"); ax.set_ylabel("$I_\\nu(x)$")
    ax.set_title("Modified Bessel functions, first kind")
    ax.legend(); _save(fig, "besseli.png")


def plot_besselk():
    x = np.linspace(0.05, 5, 600)
    fig, ax = plt.subplots()
    ax.plot(x, sp.k0(x), label="$K_0(x)$")
    ax.plot(x, sp.k1(x), label="$K_1(x)$")
    ax.set_ylim(0, 3)
    ax.set_xlabel("x"); ax.set_ylabel("$K_\\nu(x)$")
    ax.set_title("Modified Bessel functions, second kind")
    ax.legend(); _save(fig, "besselk.png")


def plot_hankel():
    x = np.linspace(0.2, 20, 800)
    h = sp.hankel1(0, x)
    fig, ax = plt.subplots()
    ax.plot(x, h.real, label=r"$\mathrm{Re}\,H_0^{(1)}(x)$")
    ax.plot(x, h.imag, label=r"$\mathrm{Im}\,H_0^{(1)}(x)$")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlabel("x"); ax.set_ylabel("$H_0^{(1)}(x)$")
    ax.set_title("Hankel function of the first kind")
    ax.legend(); _save(fig, "hankel.png")


def plot_gamma():
    x = np.linspace(-4.5, 5, 2000)
    x = x[np.abs(x - np.round(x)) > 2e-2]  # drop points near the integer poles
    g = sp.gamma(x)
    fig, ax = plt.subplots()
    ax.plot(x, g, ".", ms=1.5)
    ax.set_ylim(-20, 20); ax.axhline(0, color="gray", lw=0.6)
    for p in range(0, -5, -1):
        ax.axvline(p, color="gray", lw=0.5, ls=":", alpha=0.5)
    ax.set_xlabel("x"); ax.set_ylabel(r"$\Gamma(x)$")
    ax.set_title("Gamma function")
    _save(fig, "gamma.png")


# ---------------------------------------------------------------------------
# Tier B — true accuracy of fortran-bessels vs an mpmath reference
# ---------------------------------------------------------------------------

def _have_fpm():
    try:
        subprocess.run(["fpm", "--version"], cwd=ROOT,
                       capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def dump(name, order=None):
    """Run the Fortran helper and return (x, f) arrays of the library's output."""
    cmd = ["fpm", "run", "--example", "dump_values", "--profile", "release",
           "--", name]
    if order is not None:
        cmd.append(str(order))
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True)
    xs, fs = [], []
    for line in out.stdout.splitlines():
        line = line.strip()
        if not line or line[0] not in "+-.0123456789":
            continue  # skip any stray fpm/runtime chatter
        try:
            a, b = line.split(",")
            xs.append(float(a)); fs.append(float(b))
        except ValueError:
            continue
    (DATA / f"{name}.csv").write_text(out.stdout)
    return np.asarray(xs), np.asarray(fs)


def ulp_error(f, ref, envelope=False):
    """Error in ULPs of real64. For monotonic functions this is the ordinary
    relative error |f-ref|/|ref|. For oscillatory functions (envelope=True) the
    denominator is the local oscillation amplitude (a rolling max of |ref|), so
    the metric measures error relative to the wave amplitude instead of blowing
    up at the function's zeros."""
    from scipy.ndimage import maximum_filter1d
    ref = np.asarray([float(r) for r in ref])
    denom = maximum_filter1d(np.abs(ref), size=51, mode="nearest") if envelope \
        else np.abs(ref)
    denom = np.maximum(denom, np.finfo(float).tiny)
    err = np.abs(f - ref) / denom / EPS
    return np.where(err == 0.0, np.nan, err)  # drop exact matches on the log axis


def accuracy_plot(fname, title, members, envelope=False):
    """members: list of (label, fortran_name, mpmath_callable). One curve each."""
    import mpmath as mp
    mp.mp.dps = 40
    fig, ax = plt.subplots()
    for label, name, fn in members:
        x, f = dump(name)
        ref = [fn(mp.mpf(float(xi))) for xi in x]
        ax.semilogy(x, ulp_error(f, ref, envelope), ".", ms=2.5, label=label)
    ax.axhline(1.0, color="r", lw=0.9, label="1 ULP")
    ax.set_xlabel("x")
    ax.set_ylabel("error / amplitude  [ULP]" if envelope
                  else "relative error  [ULP]")
    ax.set_title(title); ax.legend()
    _save(fig, fname)


def tier_b():
    import mpmath as mp
    accuracy_plot("besselj_accuracy.png", "First kind J — accuracy vs mpmath", [
        ("$J_0$", "besselj0", lambda x: mp.besselj(0, x)),
        ("$J_1$", "besselj1", lambda x: mp.besselj(1, x)),
    ], envelope=True)
    accuracy_plot("bessely_accuracy.png", "Second kind Y — accuracy vs mpmath", [
        ("$Y_0$", "bessely0", lambda x: mp.bessely(0, x)),
        ("$Y_1$", "bessely1", lambda x: mp.bessely(1, x)),
    ], envelope=True)
    accuracy_plot("besseli_accuracy.png", "Modified I — accuracy vs mpmath", [
        ("$I_0$", "besseli0", lambda x: mp.besseli(0, x)),
        ("$I_1$", "besseli1", lambda x: mp.besseli(1, x)),
    ])
    accuracy_plot("besselk_accuracy.png", "Modified K — accuracy vs mpmath", [
        ("$K_0$", "besselk0", lambda x: mp.besselk(0, x)),
        ("$K_1$", "besselk1", lambda x: mp.besselk(1, x)),
    ])
    accuracy_plot("gamma_accuracy.png", "Gamma — accuracy vs mpmath", [
        (r"$\Gamma$", "gamma", lambda x: mp.gamma(x)),
    ])


def _save(fig, fname):
    fig.tight_layout()
    fig.savefig(MEDIA / fname, bbox_inches="tight")
    plt.close(fig)
    print(f"  {fname}")


if __name__ == "__main__":
    print("Generating fortran-bessels documentation plots...")
    print(" curves:")
    plot_besselj(); plot_besseljn(); plot_bessely()
    plot_besseli(); plot_besselk(); plot_hankel(); plot_gamma()

    if _have_fpm():
        try:
            import mpmath  # noqa: F401
        except ImportError:
            print(" accuracy: SKIPPED (pip install mpmath)")
        else:
            print(" accuracy:")
            tier_b()
    else:
        print(" accuracy: SKIPPED (fpm toolchain not found)")

    print(f"\nDone. Plots written to {MEDIA}/")
