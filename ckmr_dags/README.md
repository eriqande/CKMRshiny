# CKMR DAGs

Start with `ckmr_dags.qmd`, the growing document for DAGs across course examples.
Its first example describes one candidate mother-offspring pair in the
female-only Acme model: observed ages/years, deterministic birth-year and
age-at-birth bookkeeping, focal and total expected reproductive outputs, true
kinship, and detected kinship. The two reproductive-output nodes stay visible.

The four top nodes are constrained to the same rank and the order
`A_Ju, Y_Ju, A_Ad, Y_Ad`. Invisible DOT edges impose that layout; they are not
model edges. Sex is conditioned on being female and is not another node.

## Build

From the repository root, with the base environment active:

```sh
python3 ckmr_dags/build_ckmr_dags.py --document
```

Requirements: Python 3, Graphviz (`dot`), `dot2tex`, a LaTeX installation with
TikZ/standalone/Latin Modern, Quarto, and R with knitr/rmarkdown. All executables
must be on PATH. The build script itself uses only Python's standard library.
The R example uses base R and knitr, and its numerical assertions execute when
Quarto renders the walkthrough.

The pipeline is deliberately external to Quarto:

1. `acme.dot` defines nodes, edges, LaTeX labels (`texlbl`), colors, and ranks.
2. Graphviz selects the layout. The build checks the top-row constraints.
3. dot2tex generates TikZ, and pdflatex typesets actual LaTeX node labels.
4. Quarto includes the resulting vector PDF as an ordinary external figure.

There is no fenced `dot` block. The generated TikZ is not edited by hand.

Outputs:

- `output/pdf/ckmr_dags.pdf`: walkthrough, factorization, marginalization, and R code.
- `output/pdf/acme-dag.pdf`: standalone graph for inspection or reuse.

`figures/acme-dag.pdf` is the generated figure included by the document. Build
logs, generated TeX, and Graphviz layout JSON go to `../tmp/pdfs/acme-dag`.
All generated files are ignored by Git; retain the DOT, QMD, and build script.
Rebuild only the graph with `python3 ckmr_dags/build_ckmr_dags.py`.

## Conventions and assumptions

- Light blue fill and double outlines: pair-member covariates conditioned upon.
- Gray fill: other observed quantities. White: unobserved. This does not imply conditional
  randomness: both reproductive-output nodes are deterministic in this model.
- Dark pink fill: estimated parameters `N_{0,f}` and `lambda`.
- Pink outline: a value or conditional law depends on the parameter vector.
  In particular, the observed value of `K*` stays fixed even though its
  conditional probability depends on parameters.
- Dashed outline/incoming arrows: data-only bookkeeping intermediates.
- Light blue diamonds: supplied constants or modeling choices (`y0`, maturity age,
  false-positive rate `alpha`, and false-negative rate `beta`).

The aim is to keep alpha effectively zero. Mark's Acme example takes both
alpha and beta to be zero. Lambda here is the instantaneous exponential growth
rate (annual multiplier `exp(lambda)`). Abundance means mature female abundance.
Lethal sampling follows breeding, so capture in the offspring's birth year is
allowed. `Ju` is a candidate-offspring role, not necessarily a juvenile at capture.

The model agrees with `generic_lglk_MOP_ideal_mammal` in the course repository's
`r-scripts/ckmr_funs.r`, as used by `fit_delfi_A.r`. The existing
`../tidy_rmd_examples/fit_delfi_A.Rmd` provides the data-preparation walkthrough.
Pairs are represented through individual marginals; the graph does not assert
independence across overlapping pairs. Grouping comparisons and applying Mark's
Poisson pseudolikelihood is described separately from the pair's Bernoulli model.
