# Validation report: opt-in BDI extension

Validation session: 2026-09-17.
Base CAFE5 commit: `b9e3b2e55a2412fca94f69ae42a2028cfeca3e95`.
Environment: x86-64 WSL2, GCC 15.2, CMake 4.2.3, C++11, OpenMP.

## Implemented and checked

- Analytical equal-duplication/loss plus innovation transition kernel, including
  exact zero-rate boundaries and positive transitions out of zero.
- Joint global lambda/nu optimization, deterministic multistart, and zero-rate
  boundary fits; fixed rates are also supported.
- Explicit zero-inclusive root law, marginal root likelihood, scaled pruning,
  repeated-pattern reuse, and observed-family ascertainment correction.
- Nodewise marginal ancestral posteriors, including root zero.
- Unbounded branch simulation under the same model and sampling condition.
- Opt-in experimental fixed-parameter family Monte Carlo checks.
- Root-distribution sensitivity and both evaluation and refitting checks for
  finite count truncation.

This mode uses CAFE5's tree representation/parser and optimizer, with a separate
likelihood/reconstruction/simulation path. It does not retrofit innovation into
the legacy gamma/error-model machinery. Gamma mixtures, jointly fitted error
models, branch-specific innovation, separate duplication/loss rates, a unique
family-origin process, and calibrated branch significance are not implemented.

## Numerical checks

The independent Python suite compared the C++ transition kernel to exponentiation
of an independently constructed CTMC generator, with a larger state space for the
reference calculation. Maximum entry errors in the tested cases were below
`2.2e-14`, including lambda*t above the legacy saturation regime. Zero rates,
zero time, near-zero lambda, row mass, moments and the semigroup identity passed.
A high-immigration Poisson case also checks that initial tail underflow does not
zero out the distribution near its mode.

On `((A:0.3,B:0.3):0.4,C:0.7)`, independent generator-based pruning agreed with
conditional marginal likelihoods to within `6.3e-15`. Root and internal-node
marginal posteriors also agreed. The total optimized-engine objective was checked
against the sum of separately computed family likelihoods.

An independent **base R uniformization** calculation agreed to `4.6e-15` on low
states. This provides a second reference implementation without calling the
production BDI recurrence. See `scripts/innovation/check_kernel.R`.

The existing CAFE suite passed **206 tests and 582 assertions**, with no failures.
The build uses warning-specific exceptions for unchanged upstream code promoted
to errors by GCC 15; warnings are not globally disabled.

Machine-readable checks: [numerical_checks.json](validation/numerical_checks.json).

## Simulation and estimation

Twelve independently seeded datasets, each with 2,500 observed families on a
six-tip tree, were generated at lambda=0.2, nu=0.3, with a Poisson(1) root law.
Each dataset was fitted from multiple starts with both rates free and with exact
zero-rate boundary alternatives. Mean recovered estimates:

| Parameter | Generating value | Mean fitted value |
|---|---:|---:|
| lambda | 0.2 | 0.19990698 |
| nu | 0.3 | 0.30270671 |

All final truncation checks passed. This is recovery in one generating regime,
not a comprehensive identifiability or coverage study. Individual estimates are
in [recovery.tsv](validation/recovery.tsv).

Additional simulation checks verified moments, reproducibility within this C++
runtime, and 3,000 observed families generated from an exactly zero root.

For 1,500 independent null families scored against 5,000 fixed-parameter null
simulations, the fraction with Monte Carlo p<=0.05 was **0.03533**. The discrete
likelihood score and finite reference sample can make this conservative. This
limited experiment does not establish calibration after estimating parameters
from the tested data. Accordingly, the option remains explicitly experimental;
no branch p-values or biological significance calls are reported.

## Primary wasp test: Vespidae_dated_primary + N11

The selected data comprise **13,836 HOGs, 15 species, 2,739 distinct tip-count
patterns**, and a maximum observed count of **145**. Counts are numbers of gene
IDs in HOG membership cells, not expression counts. No TE or large-family
filter was applied. All 13,836 rows are observed in at least one selected taxon.

The tree's rooted topology matches the source OrthoFinder N11 clade. In particular,
Ancistrocerus is included in both. N13 excludes Ancistrocerus and would introduce
structural out-of-scope zeros if used on this tree. N10 matches the ancestor needed
for the initially requested tree with the two external outgroups; N12 is its
non-Vespidae daughter. Node labels are specific to this source OrthoFinder run.

The primary tree is the actual `Vespidae_dated_primary.nwk` file, not a pruning of
the separately dated `Vespidae_with_outgroups_dated_primary` file.
Source checkout: `tomwenseleers/evodevo_waspcastes`, local HEAD
`2cbde4cee15502dd31df434f6d38009073f06689`. Exact file hashes, rather than only HEAD,
identify the analyzed inputs:

| Input | SHA256 |
|---|---|
| Source N11.tsv | `ab9a23a3be518c0c6d45146c3eed3891d6aba62e24b7cb3ed62d1844b761d49f` |
| Source Vespidae_dated_primary.nwk | `d8f3dd6387f3222bffd2191ffae3c936f727e17760026370083c5ec69caa55f3` |
| Generated N11 CAFE count table | `c594ea94e7aa7741602a39dc16b502187a1efca42a174fc543a6efb46c019562` |

Full input audit: [input_provenance.json](validation/input_provenance.json).

With a fixed **Poisson(1) root distribution**, joint estimates are

- lambda = **0.00140003266** per copy per million years;
- nu = **0.00121697482** per family per million years;
- conditional negative log likelihood = **114467.3778354**.

The supplied dating calibration and QC tables express branch lengths in Ma, so
these fitted rates are per million years.

Both positive-rate starts converged to essentially the same estimates. For the
same root law and likelihood formulation, the best nu=0 boundary had negative
log likelihood **136021.7334040**, and the best lambda=0 boundary had
**189463.6062871**. These are likelihood comparisons, not calibrated significance
results. The nu=0 process reduces to ordinary birth–death, but its marginal,
ascertainment-corrected objective is not the legacy CAFE maximized-root objective;
one should not compare their reported scores as though they were identical.

The model assigned a **marginal MAP root count of zero to 4,461 families**; 3,504
had posterior probability of root zero >0.95. These posterior results are
conditional on the fitted model and specified root law and do not independently
establish family origin histories.

## Sensitivity

| Fixed Poisson root mean | lambda | nu | Negative log likelihood |
|---:|---:|---:|---:|
| 0 (all families start absent) | 0.000811070 | 0.004175069 | 126754.256516 |
| 0.5 | 0.001456248 | 0.000866346 | 111535.429444 |
| 1 | 0.001400033 | 0.001216975 | 114467.377835 |
| 2 | 0.001438671 | 0.001396933 | 123325.959302 |

The zero-root row is an extreme sensitivity model, not a preferred biological
assumption. Innovation estimates depend materially on the root law. The mean-one
fit is a reproducible reference, not a data-selected optimum for the root law.
The lower score at mean 0.5 indicates that Poisson(1) should not simply be assumed
adequate. More flexible root laws, observation-error sensitivity and biological
model adequacy deserve further investigation before interpreting a rate as a
robust biological estimate.

For every listed fit, raising the cap from **180 to 360** changed total negative
log likelihood by zero at the reported double precision. Independently refitting
the mean-one model at **K=240** recovered the same printed lambda, nu and likelihood;
its final **K=480** check also passed. This supports truncation convergence for
these fits. The complete sensitivity table is
[empirical_sensitivity.tsv](validation/empirical_sensitivity.tsv).

## Reproduce

```bash
python3 scripts/innovation/prepare_wasps.py /path/to/evodevo_waspcastes validation/wasps
./build/cafe5 --innovation \
  -t validation/wasps/Vespidae_dated_primary.nwk \
  -i validation/wasps/N11_Vespidae_counts.tsv \
  --root-mean 1 --max-count 180 --threads 3 --starts 2 \
  -o validation/wasps/N11_root1
```

Repeat with root means 0, 0.5 and 2; repeat mean 1 with `--max-count 240`.
Generated empirical inputs, full ancestral tables, simulations and executables are
kept locally under ignored `validation/` and `bin/` directories. The repository
contains the reproducible scripts and compact audit/validation summaries.
