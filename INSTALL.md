# Building and installing `bdpaft`

## 1. Once-only setup

```r
install.packages(c("devtools", "roxygen2", "testthat", "Rcpp", "RcppArmadillo"))
```

## 2. First-time package generation

After extracting this scaffold into your existing `bdpaft/` repo:

```r
setwd("path/to/bdpaft")
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
# R CMD INSTALL bdpaft
```

First compile takes ~30-60 s (RcppArmadillo headers are large).

## 4. Run tests

```r
devtools::test()                # quick smoke tests (<1 min)
```

## 5. Full CRAN check (before submission, much later)

```bash
R CMD build bdpaft
R CMD check --as-cran bdpaft_0.1.0.tar.gz
```

Expect a few notes the first time around — they will need cleanup before CRAN
submission but are fine for an internal/lab release.

## 6. Edit the placeholder author info

Before pushing, update three placeholders:

- `DESCRIPTION`: `Authors@R` — set your real name and email
- `LICENSE`: `COPYRIGHT HOLDER` — same
- `LICENSE.md`: copyright line — same

## 7. Sanity-check the C++ source

The bundled `src/bdpaft.cpp` was derived from your standalone
`bdpaft.cpp`. If your production version is actually `bdpaft_v2.cpp`,
diff the two and port any v2-only changes into `src/bdpaft.cpp` before
`devtools::install()`:

```bash
diff -u /path/to/bdpaft.cpp /path/to/bdpaft_v2.cpp
```

Any differences should be applied to `src/bdpaft.cpp`. Watch in particular
for the function signature: in the package version the exported function is
named `bdpaft_cpp` (renamed from `bdpaft` to avoid colliding with the R
wrapper).

## 8. Optional: keep your analysis repo separate

This package contains only the algorithm. Paper-specific analysis scripts
(R00..R04, output/) should live in a sibling repo such as `bdpaft-cathgen`.

## What this scaffold does NOT include yet

- Vignette (will be needed before Biostatistics submission for code-review)
- Plotting methods
- Example data (the smoke tests synthesize data on the fly)
- `pkgdown` site

These can all be added incrementally without breaking what is here.
