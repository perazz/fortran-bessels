# Documentation

The doc site is built with [Doxygen](https://www.doxygen.nl) (1.13.2) and the
[doxygen-awesome](https://github.com/jothepro/doxygen-awesome-css) theme, and is
deployed to GitHub Pages by `.github/workflows/deploy-docs.yml` on every push to
`main`. The live site is at <https://perazz.github.io/fortran-bessels>.

Layout:

```
project/doxygen/Doxyfile        Doxygen configuration
project/doxygen/header.html     HTML header (dark-mode toggle)
project/doxygen/doxygen-awesome*.{css,js}   vendored theme (MIT, jothepro)
doc/mainpage.md                 doc-site landing page
doc/{besselj,bessely,besseli,besselk,hankel,gamma}.md   per-family pages
doc/generate_plots.py           regenerates the figures in doc/media/
doc/media/*.png                 committed figures (curves + accuracy)
```

## Build the site locally

Prerequisites: `doxygen` (>= 1.9) and `graphviz` (for the `dot` graphs).

```bash
brew install doxygen graphviz          # macOS
cd project/doxygen && doxygen
open html/index.html
```

The committed PNGs in `doc/media/` mean a plain `doxygen` run needs no Python.

## Regenerate the plots

Prerequisites: Python 3 with `numpy`, `scipy`, `matplotlib`, `mpmath`.

```bash
pip install numpy scipy matplotlib mpmath
python doc/generate_plots.py           # writes doc/media/*.png
```

Re-run this and commit the updated PNGs only when the figures need to change.

### How the accuracy plots are produced

The `*_accuracy.png` figures compare this library's own output against an
arbitrary-precision `mpmath` reference. `generate_plots.py` drives the small
Fortran helper `example/dump_values.f90`, which prints `x,f(x)` for a named
function:

```bash
fpm run --example dump_values -- besselj0
```

`generate_plots.py` captures that CSV (into `doc/media/data/`, gitignored),
evaluates the same points with `mpmath` at high precision, and plots the relative
error in ULPs of `real64`. Building the helper requires an `fpm`/`gfortran`
toolchain; if it is unavailable, the accuracy plots are skipped and only the
function-shape curves are produced.
