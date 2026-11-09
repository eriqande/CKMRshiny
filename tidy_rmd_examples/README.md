# Tidy CKMR notebooks

Start with [fit_delfi_A.Rmd](fit_delfi_A.Rmd), a static R Markdown companion to
Mark's Acme Archipelago example from the 2026 Hawaii course. Potential
offspring are selected by birth year (2000 onward by default), including
individuals sampled as adults. Pair entries are individual IDs. Preparation is shown inline, with
inspectable tables, figures, and mathematics. Preparation stays inline;
likelihood functions are introduced when repeated evaluation becomes useful.

Continue with [fit_delfi_S.Rmd](fit_delfi_S.Rmd) for **Sanctuary Sound**. This
notebook currently covers the supplied growth/length model, explicit length bins,
two-epoch age-length calibration, and a tidy table of observed comparison
categories. A worked category expands both uncertain ages into 400 terms and
shows their contributions in heat maps, including post-capture survival. It does
not yet fit the demographic model or implement the efficient full likelihood.
Known ages are used for calibration and omitted from the kin-analysis sample.
Exact unordered-pair counts and Mark's same-cell approximation are kept separate.

Sanctuary's example-specific reusable functions live together in
[R/tidy-delfi-S.R](../R/tidy-delfi-S.R), in order of use. Its YAML parameters
default to `../.hawaii-ckmr-26/datasets/samp_delfi_S.rda` and
`../R/tidy-delfi-S.R`. As with Acme, paths resolve from the Rmd folder. The
optional collapsed checks compare calibrated probabilities and Mark-compatible
category counts with the original functions, and verify the one-row mixture.
Set `check_mark_code: false` to render without those checks. Small analytic helper
checks can be run from the project root with `Rscript tests/check-tidy-delfi-S.R`.

Open the Rmd in RStudio and run its chunks in order or click **Knit**. Required
packages are `dplyr` (1.1.0 or newer), `tidyr`, `ggplot2`, `tibble`, `rmarkdown`,
and `knitr`. Rendering HTML also requires Pandoc, included with RStudio.

Change `data_path` in the first data-loading chunk, or edit the `data_path` YAML
parameter. The current default is `../.hawaii-ckmr-26/datasets/samp_delfi_A.rda`, relative
to the Rmd's folder. Absolute paths work too. The simulated data are not
copied into this folder; participants need their own authorized copy of Mark's
course repository.

To render from the CKMRshiny project root:

```r
rmarkdown::render(
  "tidy_rmd_examples/fit_delfi_A.Rmd",
  output_dir = "tidy_rmd_examples/output"
)
```

Here `output_dir` is relative to the console's working directory (the project
root), while relative data paths inside chunks resolve from the Rmd's folder.
To override the data location explicitly:

```r
rmarkdown::render(
  "tidy_rmd_examples/fit_delfi_A.Rmd",
  params = list(data_path = "/absolute/path/to/hawaii-ckmr-26/datasets/samp_delfi_A.rda",
                check_mark_code = FALSE),
  output_dir = "tidy_rmd_examples/output"
)
```

The original-code cross-check compares preparation and likelihood values with
Mark's original functions, and checks the fitted parameter estimates. It uses
`mvbutils`, `offarray`, `atease`, and `kinsimmer2` only when available and when
`check_mark_code` is true. It sources `r-scripts/ckmr_funs.r` from the repository
containing the chosen dataset; it does not run Mark's interactive driver script.

The final optional section introduces **RTMB**, using an ordinary parameter
vector and explicit data inputs. Set `use_rtmb: false` in the YAML header to
skip execution, or leave it true to run when RTMB is installed. The document
does not install RTMB automatically. It explains fixed data versus taped
arithmetic, structural zeros, supported probability functions, parameter-dependent
checks, and rebuilding the tape after data changes. It checks objective values,
gradients, and fitted parameters against the ordinary R version.

Uncertainty calculations first reuse the manual delta method with an AD Hessian,
then reproduce it with `ADREPORT()` and `sdreport()`, explicitly supplying the
same Hessian. Both log abundance and abundance are reported; confidence limits
are constructed on the log scale and exponentiated, as in the earlier analysis.
The section also explains named parameter lists and `getAll()` as alternative
interfaces without requiring them for this example. It was verified with RTMB 2.0.

Generated HTML files are kept in `output/`, which is ignored by Git.
