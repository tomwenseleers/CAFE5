# Gamma variation, observation error and branch-tail tests

This extends the opt-in innovation model. Gamma variation affects **duplication
and loss only**, as selected for this analysis. The innovation rate remains one
global nu. Both fixed and estimated annotation error are supported.

## Generative model

For family f, draw one rate multiplier G_f from a discrete approximation to
Gamma(shape=alpha, rate=alpha), normalized to unit mean. This single family category
applies to **every branch** of its tree:

\[
q^{(f)}_{n,n+1}=\lambda G_f n+\nu,\qquad
q^{(f)}_{n,n-1}=\lambda G_f n.
\]

The implementation uses CAFE5's existing median-category discretization, not an
independent category draw on each branch. Smaller alpha produces greater rate
variation. The supported numerical shape interval is [0.05,100]; a fit near a
bound is flagged. Compare different category counts and the homogeneous model
when evaluating model adequacy. Gamma variation does not remove the mean upward
drift of the equal-duplication/loss innovation process.

True tip count N_s is latent. The observed count Y_s follows an error law
E(y|n). The global symmetric one-copy model is

\[
E(n-1|n)=\epsilon,\quad E(n|n)=1-2\epsilon,\quad
E(n+1|n)=\epsilon\quad(n\ge1),
\]

and `E(0|0)=1-epsilon`, `E(1|0)=epsilon`. Thus epsilon is the probability of
**each** one-copy error at positive counts, not the total positive-count error
probability (which is 2*epsilon). Valid epsilon values are [0,0.5). False positives
at true zero are explicitly represented. Error is independent among tips given
true counts. This is an observation-error model combined with the gamma mixture;
it does not introduce separate latent error-rate categories or species-specific
error parameters.

Fixed CAFE-style error files are also supported:

```text
maxcnt: 200
cntdiff: -1 0 1
0 0 0.95 0.05
1 0.05 0.90 0.05
```

Rows must sum to one, deviations must be distinct integers, true count zero must
be specified, and negative observed counts must have zero probability. Omitted
rows use the preceding distribution; the last row also extends beyond `maxcnt`
for unbounded simulation and larger-state truncation checks. The header is a
compatibility label, not a cap on possible true or observed counts. Such a fixed
file cannot be combined with fixed or estimated epsilon.

## Correct observed-family likelihood

Let L_fk be the family likelihood in rate category k, integrating root states and
all latent ancestral/tip counts. At each tip the pruning vector is E(y_s|n), not
a point mass at the observed count. Let Z_k denote the likelihood for **all observed
counts zero**, including the error law. With pre-selection category weights w_k,

\[
\Pr(Y_f=y_f\mid Y_f\ne0)
=\frac{\sum_k w_kL_{fk}}{1-\sum_k w_k Z_k}.
\]

It is generally incorrect to average separately conditioned category likelihoods
with the original weights. Equivalently, observed families have category weights
proportional to `w_k*(1-Z_k)`. Simulation draws a category and latent history,
then applies observation error, then rejects an all-zero observed vector. It
redraws the category as well as the history after rejection.

Category posterior probabilities given a nonzero observation are proportional
to `w_k*L_fk`. Ancestral posteriors average category-specific posteriors with these
weights. The same procedure estimates latent **true** counts at tips, so terminal
branch changes do not equate an observed annotation error with an evolutionary
change.

## Fit and simulation commands

```bash
./bin/cafe5 --innovation \
  -t validation/wasps/Vespidae_dated_primary.nwk \
  -i validation/wasps/N11_Vespidae_counts.tsv \
  --root-mean 1 --gamma-cats 3 --estimate-epsilon \
  --max-count 180 --iterations 1000 --starts 2 --threads 6 \
  -o validation/mixtures/N11_gamma_error
```

`--gamma-cats` defaults to 1. When it exceeds 1, alpha is estimated unless
`--alpha` fixes it. Error defaults to zero; use exactly one of `--epsilon E`,
`--estimate-epsilon`, or `--error-model FILE`. Joint estimation of lambda, nu,
alpha and epsilon is supported. Rates use log coordinates and epsilon uses a
logit transform into (0,0.5). Exact zero-lambda, zero-nu and zero-epsilon faces are
also fitted. The root law is still supplied explicitly and fixed.

Simulation requires fixed values for every active model parameter. It uses the
same discrete gamma approximation and observation law as inference.

The result table adds alpha, epsilon, category count, error-file path and an
alpha-bound indicator. `_categories.tsv` records category weights and lambda
multipliers. `_families.tsv` also records posterior category probabilities.
`_branch_statistics.tsv` records posterior mean change for each
original family ID and each non-root node. All previous likelihood and
truncation checks remain active. A nonconverged or truncation-failing fit returns
nonzero status and must not be silently accepted in a bootstrap.

## What the branch tests mean

For family f and parent-to-child branch b, the statistic is

\[
T_{fb}(Y_f)=E[N_{child}-N_{parent}\mid Y_f;\theta].
\]

It is calculated from marginalized **true-count** posteriors. The null population
is an exchangeable observed family generated by the global BDI+gamma+error model
on the same tree, using the same root law and inclusion rule. Upper and lower tails
assess unusually large and unusually small inferred changes on that branch.
This is **not** the null hypothesis that there was no copy-number change, and it
is not a branch-specific rate-ratio likelihood-ratio test. It differs from CAFE's
legacy Viterbi branch quantities. A significant model-tail result does not
establish selection, molecular innovation, or a unique family origin.

These are population-level family/branch checks, not tests conditioned on a known
ancestral family size or a known family-specific rate. Consequently, an unusual
family-wide history or a misspecified root law may yield extreme branch scores.
Inspect whole-family fit and root-law sensitivity when interpreting the results.

### Known/fixed-parameter Monte Carlo mode

With genuinely specified parameters, simulate B independent observed families,
reconstruct them with those same parameters, and compare each observed branch
score against the corresponding simulated branch scores:

\[
p^+_{fb}=\frac{1+\sum_{r=1}^{B}I(T_{rb}^*\ge T_{fb})}{B+1},\quad
p^-_{fb}=\frac{1+\sum_{r=1}^{B}I(T_{rb}^*\le T_{fb})}{B+1}.
\]

The two-sided value is `min(1,2*min(p_upper,p_lower))`. Under the specified null,
observed and simulated statistics are exchangeable. Random ranks would be
uniform; including ties in the tail makes the reported ranks conservative. Thus
each one-sided p-value is super-uniform. The union bound proves conservativeness
of the doubled-minimum two-sided value. This argument requires identical scoring,
conditioning and fixed parameters; it does not justify treating fitted parameters
as known without qualification.

```bash
python3 scripts/innovation/branch_bootstrap.py ./bin/cafe5 \
  -t tree.nwk -i counts.tsv -o branch_fixed \
  --root-mean 1 --gamma-cats 3 --alpha 1.3 --epsilon 0.08 \
  --lambda 0.2 --nu 0.3 --fixed-parameters --replicates 9999
```

### Whole-dataset bootstrap with refitting (default)

1. Fit the supplied data under its requested parameter constraints.
2. Simulate a dataset with the same number of observed families at the fitted
   parameters, including root/category draws, evolution and annotation error.
3. Refit all parameters that were free in the original fit; keep originally fixed
   parameters fixed. Reconstruct all branch statistics.
4. Compare each observed family index with the corresponding exchangeable family
   index in each simulated dataset. This supplies B replicate statistics per
   hypothesis, not an incorrectly inflated B*number-of-families replicate count.
5. Repeat and calculate the same finite-replicate tail estimates.

```bash
python3 scripts/innovation/branch_bootstrap.py ./bin/cafe5 \
  -t validation/wasps/Vespidae_dated_primary.nwk \
  -i validation/wasps/N11_Vespidae_counts.tsv -o branch_refitted \
  --root-mean 1 --gamma-cats 3 --estimate-epsilon \
  --max-count 180 --iterations 1200 --starts 2 --threads 6 \
  --replicates 999 --seed 7142
```

This is a **plug-in parametric bootstrap**, so its composite-null calibration is
approximate, not guaranteed by the fixed-parameter rank proof. The full refitting
pipeline must be checked using independent outer simulations. Any failed fit,
truncation failure, or simulated count exceeding the cap aborts the run; no
replicate is silently dropped. Rerun with suitable computational settings.

A refitted bootstrap requires B complete model fits and is considerably more
expensive than fixed-parameter simulation. The implementation does not silently
switch to the faster method.

## Multiple tests and Monte Carlo precision

`branch_tests.tsv` includes both tails, two-sided p-values, Benjamini–Yekutieli
FDR adjustment and Holm FWER adjustment across **all supplied family–branch
hypotheses**. BY and Holm accommodate arbitrary dependence when the marginal
p-values are valid; they cannot correct miscalibration of the underlying model or
an approximate bootstrap. No prefilter based on observed significance is applied.

The minimum two-sided p-value is `2/(B+1)`. Small B may make adjusted significance
impossible across thousands of hypotheses. A 39-replicate run is a coarse pipeline
calibration test, not a discovery analysis. Even 999 replicates may be inadequate
for stringent multiplicity control over the complete N11 matrix. Plan B from the
required tail resolution and computational budget; do not report zero p-values
or infer extreme significance from an inadequate simulation count.

`method.json` records the exact null, mode, resolution and model parameters.
`bootstrap_fits.tsv` records every refitted parameter vector and diagnostic;
replicate datasets, reconstructions and logs are retained for audit.

## Validation and scientific limits

Run `validate_mixtures.py` for an independent generator-matrix likelihood,
observation-error equivalence, simulation frequencies and joint parameter recovery.
Run `validate_branches.py` for fixed-parameter and fully refitted outer calibration
experiments. The numerical references do not reuse the production BDI recurrence.
See `innovation_mixture_validation.md` for measured results and their limits.

Innovation, annotation error, gamma variation and the root distribution can
partly compensate for one another. A converged fit alone does not establish
identifiability. Fixed error models derived from external annotation validation
remain useful. Confidence intervals, broader null/power studies, category-count
sensitivity, and empirical model-adequacy checks are additional analyses; the
presence of a bootstrap implementation does not establish universal calibration.
