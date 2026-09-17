# Gamma/error and branch-test validation, 2026-09-18

The extension was developed and tested in WSL2 on the `innovation-bdi` branch.
It retains the global innovation parameter and applies discrete gamma variation
only to duplication/loss. See [model and test definitions](innovation_mixtures_and_branches.md).

## Numerical and simulation checks

An independent three-category gamma+error calculation used SciPy generator-matrix
exponentials, explicit observation-emission vectors and a zero-inclusive Poisson
root law. Conditional log likelihoods agreed with the C++ implementation to
**2.7e-15**. Root and category posterior probabilities also agreed. The existing
CAFE gamma-quantile approximation agreed with normalized SciPy median quantiles
to **4.9e-9** in the tested case.

A fixed CAFE-format error file and the equivalent epsilon model produced the
same likelihood. An independent 40,000-family simulation checked selected
observed-pattern frequencies against their explicit mixture probabilities.
The lambda=0 boundary reduces to the same pure-immigration likelihood regardless
of gamma category count; alpha is then unidentifiable and is not needlessly
estimated. Fixed user-specified alpha values remain recorded.

Three joint-recovery experiments used 4,000 observed families each, with lambda,
nu, alpha and epsilon all fitted simultaneously:

| Parameter | Generating value | Mean recovered |
|---|---:|---:|
| lambda | 0.2 | 0.196400 |
| nu | 0.3 | 0.301562 |
| alpha | 1.3 | 1.402777 |
| epsilon | 0.08 | 0.078313 |

These are three replicates at one regime, not a comprehensive identifiability or
interval-coverage study. The unchanged homogeneous-mode numerical checks and the
**206 existing tests / 582 assertions** also passed.

Machine-readable results: [mixture_checks.json](validation/mixture_checks.json).

## Branch-tail calibration

The statistic is the posterior mean true-count change on a branch, with category,
root and observation uncertainty marginalized. The null is the global observed-
family population under the specified model; it is not 'no evolutionary change'.

**Known-parameter check:** 2,000 independent null families were tested against
3,999 simulated reference families on a four-tip tree with gamma variation and
annotation error. Rejection rates at nominal two-sided 5% were:

| Branch | Rejection fraction |
|---|---:|
| Node1 | 0.0480 |
| A | 0.0460 |
| B | 0.0480 |
| Node4 | 0.0520 |
| C | 0.0435 |
| D | 0.0440 |

The across-branch mean was **0.04692**. Finite-sample conservativeness for known
parameters follows from the exchangeable-rank argument in the model guide;
this simulation additionally checks its implementation.

**Refitted-bootstrap check:** four independently generated datasets of 300
observed families each were tested with 39 full-dataset bootstrap replicates.
Lambda, nu, alpha and epsilon were re-estimated in **all 156 replicate fits**.
Mean rejection fractions across the six branches were 0.04944, 0.04167, 0.05111
and 0.05611 for the four outer datasets, averaging **0.04958**.

This checks the actual refitting/reconstruction pipeline. It does **not** prove
calibration at all parameter values, root laws, sample sizes, branch lengths or
error distributions. There are only four independent outer datasets, and the
minimum two-sided p-value with 39 replicates is exactly 0.05. This is a limited
operating-characteristic check, not a high-resolution significance analysis.

The driver aborts on failed refits or truncation failures instead of dropping
replicates. It retains all replicate inputs, diagnostics and reconstructions.
It reports raw tails, two-sided p-values, BY FDR and Holm FWER corrections across
all requested family–branch hypotheses. The multiplicity adjustments require
valid marginal p-values; they cannot repair a misspecified model or an inadequately
calibrated plug-in bootstrap.

Machine-readable results: [branch_calibration.json](validation/branch_calibration.json).

## N11 + Vespidae_dated_primary

The same 13,836 HOGs, 15 taxa, source hashes and fixed Poisson(1) root law as the
initial report were used. Three gamma categories and an estimated global epsilon
were added. Two positive-rate starts agreed on:

| Parameter | Estimate |
|---|---:|
| Mean duplication/loss lambda, per copy per Myr | 0.05210188 |
| Global innovation nu, per family per Myr | 0.000479079 |
| Gamma shape alpha | 0.1303318 |
| Per-direction one-copy error epsilon | 0.0167541 |
| Conditional negative log likelihood | 101343.907803 |

The [multistart trace](validation/N11_gamma_error_optimization.tsv) records the
joint fit and boundary comparisons. The selected solution was re-evaluated with
the final engine to generate the linked diagnostics and posterior outputs.

The mean lambda is the mean over the pre-ascertainment gamma distribution. It is
not the rate for every family. The small alpha implies substantial heterogeneity;
its mean rate should not be interpreted as interchangeable with the homogeneous
estimate. The three lambda multipliers are approximately 0.00001098, 0.05042
and 2.94957, each with pre-selection weight 1/3.
The [fit diagnostics](validation/N11_gamma_error_results.tsv) and
[category multipliers](validation/N11_gamma_error_categories.tsv) are recorded. At positive true counts, total one-copy error probability is
2*epsilon, approximately **3.35%**.

At this solution, doubling the latent count cap from 180 to 360 changed total
negative log likelihood by zero at the reported precision. This check was also
performed independently of the multistart optimizer. The best evaluated no-error
positive-rate fit had negative log likelihood approximately 103652, and the
no-innovation fit with estimated error approximately 103289. These comparisons
are not calibrated likelihood-ratio p-values.

The estimated innovation rate differs from the earlier homogeneous, error-free
fit (0.001216975 per family per Myr). This demonstrates the practical importance
of heterogeneity and annotation error; it does not establish that the new value
is biologically definitive. Root-law and category-count sensitivity, external
annotation-error validation, and biological model adequacy remain relevant.

A full high-resolution, refitted branch bootstrap on the 13,836-family N11 dataset
has **not** been run. The implementation and small-tree calibration are available;
the empirical outputs include branch statistics, not completed branch significance
claims. B bootstrap replicates require B full empirical refits. With 15 tips there
are 28 branches and 387,408 family–branch hypotheses, so simulation resolution and
multiplicity must be planned explicitly. Do not interpret a short demonstration
bootstrap as sufficient for genome-wide discoveries.

## Reproduce the validation

```bash
python3 scripts/innovation/validate_mixtures.py ./bin/cafe5 validation/mixture_checks
python3 scripts/innovation/validate_branches.py ./bin/cafe5 validation/branch_checks
```

The branch-validation default uses four outer datasets. Continuous integration
uses one outer dataset as a smoke test, in addition to the independent likelihood,
simulation, recovery and existing CAFE tests. Neither CI success nor these finite
experiments certify universal statistical calibration.
