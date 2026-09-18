# N11 model-improvement experiment

See the [comparison and decision](comparison.md), [experiment plan](../../model_improvement_plan.md)
and [mathematical methods](../../model_improvement_methods.md).
All 13,836 N11 families are retained on the matching 15-tip tree, including large,
singleton and root-zero families. These results do not replace the archived N0
manuscript benchmark and do not provide new calibrated branch p-values.

## What was tested

- Separate duplication/loss rates; gamma still scales both and leaves innovation global.
- Separately fixed or fitted false occurrence at true count zero.
- Poisson, hurdle shifted-Poisson, hurdle shifted-NB and exact zero/one root laws.
- Individual and combined changes, gamma-resolution checks and selected second starts.
- Training on 11,166 families and scoring 2,670 held-out families, grouped by original OG.
- HOG membership, count inclusion and ultrametric marginal constraints.

## Software validation

Independent matrix-exponential checks agree with the transition kernel within
7e-15 in tested settings, including tiny unequal rates. Independent complete-tree
likelihoods agree within 8e-15. Conditional one-parameter recovery validates the
new fit parameters with other quantities known; it does not establish joint
identifiability in the empirical data. All 206 legacy tests / 582 assertions pass.

The exact observed-family simulation proposal is checked against analytic rare
immigration/error/root probabilities, multiple arrivals, independent full-tree
frequencies and gamma/error mixture frequencies. The latter also includes three
small-tree joint recovery replicates; these are not universal coverage guarantees.
The validation JSON files retain numerical results.

## Reproduction

Build with CMake as in the repository workflow. The screening, holdout preparation,
holdout scoring, report and base-R plotting scripts are in `scripts/innovation/`.
`fit_comparison.tsv` records parameters; `predictive_comparison.tsv` records full-data
predictive summaries and `per_species_predictive_checks.tsv` records tip diagnostics.
The local `validation/model_improvements/` directory additionally retains fitted
outputs, simulated count tables, commands and hash manifests. Fixed-gamma checks
must not be mistaken for refitted likelihood comparisons. Original screening runs
use their recorded frozen binary; the new conditional simulator changes random
number streams while preserving the target distribution.

No new branch bootstrap is warranted under the adequacy gate: the extended fits
remain exploratory and the existing bootstrap drivers reject these optional
parameterizations until matching refitting/reconstruction calibration is added.
