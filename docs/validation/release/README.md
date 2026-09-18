# Release evidence

- `comparison.tsv`: matched-root, matched-normalization AIC comparison.
- `gamma_resolution.tsv`: six-parameter model at verified neighboring K values.
- `benchmark_provenance.json`: exact input and engine hashes, optimization scope.
- `*_results.tsv`: score/fit diagnostics. Files named `replay` are fixed-parameter
  score verification; fitted parameter counts are given in `comparison.tsv`.
- `*_starts.tsv`, `critical_gamma_refinement.tsv`: restricted-model convergence evidence.
- `software_validation.json`: unit-test and workflow smoke-check summary.
- `extended_validation.json`: independent asymmetric-process/root/error checks.
- `rare_likelihood_validation.json`, `broad_tail_validation.json`: analytic
  checks for rare inclusion, the numerical switch, and broad count tails.
- `transition_validation.json`: independent transition-score calculation and
  known-parameter operating characteristics at one simulated regime.
- `root_validation.json`: root-law parameter recovery with other parameters fixed.

These checks establish particular numerical and software properties. They do not
certify universal identifiability, exact bootstrap size, biological model adequacy,
or a globally optimal parameter estimate. See the [mathematical specification](../../innovation_mathematics.md)
and [benchmark interpretation](../../innovation_benchmark.md).
