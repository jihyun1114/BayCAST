# Building and installing `baycast`

## 1. Once-only setup

```r
install.packages(c("devtools", "roxygen2", "testthat", "Rcpp", "RcppArmadillo"))
```

## 2. First-time package generation

After extracting this scaffold into your existing `baycast/` repo:

```r
setwd("path/to/baycast")
devtools::document()        # regenerates man/*.Rd and NAMESPACE
devtools::load_all()        # source the package without installing
```

`devtools::document()` will populate the `man/` directory from the roxygen
comments in `R/*.R`. The provided `NAMESPACE` is a hand-written placeholder
that gets overwritten on first `document()` — that is expected.

## 3. Compile & install locally

```r
devtools::install(".")
# or, equivalently, from the parent directory:
# R CMD INSTALL baycast
```

First compile takes ~30-60 s (RcppArmadillo headers are large).

## 4. Run tests

```r
devtools::test()                # quick smoke tests (<1 min)
```

## 5. Full CRAN check (before submission, much later)

```bash
R CMD build baycast
R CMD check --as-cran baycast_0.1.0.tar.gz
```

Expect a few notes the first time around — they will need cleanup before CRAN
submission but are fine for an internal/lab release.

## 6. Edit the placeholder author info

Before pushing, update three placeholders:

- `DESCRIPTION`: `Authors@R` — set your real name and email
- `LICENSE`: `COPYRIGHT HOLDER` — same
- `LICENSE.md`: copyright line — same

## 7. Sanity-check the C++ source

The bundled `src/baycast.cpp` was derived from your standalone
`baycast.cpp`. If your production version is actually `baycast_v2.cpp`,
diff the two and port any v2-only changes into `src/baycast.cpp` before
`devtools::install()`:

```bash
diff -u /path/to/baycast.cpp /path/to/baycast_v2.cpp
```

Any differences should be applied to `src/baycast.cpp`. Watch in particular
for the function signature: in the package version the exported function is
named `baycast_cpp` (renamed from `baycast` to avoid colliding with the R
wrapper).

## 8. Optional: keep your analysis repo separate

This package contains only the algorithm. Paper-specific analysis scripts
(R00..R04, output/) should live in a sibling repo such as `baycast-cathgen`.

## What this scaffold does NOT include yet

- Vignette (will be needed before Biostatistics submission for code-review)
- Plotting methods
- Example data (the smoke tests synthesize data on the fly)
- `pkgdown` site

These can all be added incrementally without breaking what is here.
