# Experimental duplication–loss–innovation model

`cafe5 --innovation` selects a separate likelihood engine for a global
copy-independent innovation rate. The default CAFE5 model is unchanged.
This is a research extension, not a validated replacement for every CAFE5
analysis. Gamma variation, fixed/estimated annotation error and model-based
branch bootstrap tests are now described in
[the extension guide](innovation_mixtures_and_branches.md). Branch-specific
rates and the legacy family/branch significance calculations remain unsupported.

## Model and interpretation

For each family independently, on every branch,

\[
q_{n,n+1}=\lambda n+\nu,\qquad q_{n,n-1}=\lambda n,
\qquad n\in\{0,1,2,\ldots\}.
\]

Rates are global, nonnegative, and expressed per tree time unit. Lambda is a
per-copy duplication rate and an equal per-copy loss rate. Nu is an acquisition
rate per family per time unit. It permits repeated acquisition and reacquisition
of the same family; it is not a unique-origin (Dollo) model and is not, by itself,
evidence for molecular de novo gene birth. HOG definitions, missed annotations,
assembly quality, horizontal transfer, and gene-family clustering affect its
biological interpretation.

Let `x=lambda*t` and `q=x/(1+x)`. One initial copy has descendant probabilities
`Pr(D=0)=q`, `Pr(D=j)=(1-q)^2*q^(j-1)` for `j>=1`. Immigration descendants have PGF

\[
H(z)=[1+\lambda t(1-z)]^{-\nu/\lambda}.
\]

Thus `I ~ NB(size=nu/lambda, prob=1/(1+lambda*t))` in R's parameterization,
and `P(i,j)` is the coefficient of `z^j` in `H(z)*G(z)^i`. The implementation
uses a positive geometric-convolution recurrence in O(K^2), rather than
repeated dense convolutions or a matrix exponential. Log probabilities initialize
the immigration row to avoid underflow before its mode. The exact boundaries are:

- `t=0`: identity;
- `nu=0`: ordinary critical birth–death, including absorbing zero;
- `lambda=0`: `N(t)=i+Poisson(nu*t)`;
- `P(0,0)=(1+lambda*t)^(-nu/lambda)` when `lambda>0`.

There is no inherited large-`lambda*t` saturation cutoff. From the generator,
`d E[N]/dt = E[(lambda*N+nu)-lambda*N] = nu`, giving
`E[N(t)|N(0)=i]=i+nu*t`. The conditional variance is
`2*i*lambda*t + nu*t + nu*lambda*t^2`. In particular this model has no
finite-mean stationary root law when `nu>0`. Its equal birth/loss constraint
imposes an upward mean trend; allowing distinct birth and death rates would be
a separate model extension.

## Root law, ascertainment and likelihood

A root distribution must be supplied explicitly. Choose `--root-mean M` for a
Poisson distribution, including zero, or `--root-prior FILE` for a finite
nonnegative distribution (two whitespace-separated columns: count and weight;
no header, unique counts; weights are normalized). Root-distribution parameters
are fixed, not jointly fitted with lambda and nu. Comparing several plausible
root laws is recommended because innovation and ancestral abundance can trade off.

Scaled pruning computes `L_f(n)=Pr(Y_f=y_f|N_root=n)`. The objective uses
**summation**, not maximization, over root counts:

\[
\ell_f=\log\sum_{n=0}^{K}\pi(n)L_f(n)
       -\log\left[1-\sum_{n=0}^{K}\pi(n)L_0(n)\right].
\]

The second term conditions on at least one observed copy among the tips. It is
computed once per rate evaluation, using `log(-expm1(log_P_zero))` for numerical
stability. Input all-zero families are rejected in this default mode. All
observed families are retained irrespective of CAFE's parsimony root filter.
`--unconditioned` explicitly removes both this correction and the corresponding
simulation rejection rule. Do not use it merely to bypass an input error.

This correction does **not** account for removing singletons, TE-associated
families, or large families. TE exclusion defines a different biological family
population and must be documented; numerical filtering by observed size requires
its own selection model. The provided empirical software tests use all observed
families in the selected HOG table, without TE or large-family filtering.

Transition entries are probabilities from the unbounded process. Finite matrix
rows are never renormalized or reflected at K. Pruning sums internal node counts
through K, so increasing K is a numerical convergence check, not an alternative
biological model. Poisson root mass omitted above K must be below 1e-8. By default,
a final evaluation with `2*K` must differ in total negative log likelihood by
less than 0.01. This checks likelihood accuracy at the fitted rates, not parameter
stability or all possible rate values. Refit with a larger K to check the latter.
`--no-truncation-check` is provided for controlled experiments and does not certify
convergence. Posterior mass at the cap is included in the reconstruction output.

Identical full patterns are weighted, and identical subtree patterns share pruning
messages. Transition matrices are reused by branch length **within one rate
evaluation only** and rebuilt whenever either rate changes. There is no cross-rate
cache that could accidentally reuse a matrix with an obsolete innovation rate.

## Build and run under WSL2

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/cafe5 --innovation --help

./build/cafe5 --innovation \
  -t Vespidae_dated_primary.nwk -i N11_Vespidae_counts.tsv \
  --root-mean 1 --max-count 180 --threads 4 --starts 3 -o vespidae
```

If CMake/SciPy are missing and sudo is unavailable, the optional
`scripts/innovation/bootstrap_wsl.sh` downloads Ubuntu packages into `.wsl-tools`
without installing system packages. Source `.wsl-tools/env.sh` afterward.
GCC, make, R (for the R check) and zlib development headers must already exist.

An omitted `--lambda` or `--nu` is estimated. Supplying it fixes that rate,
including exactly zero. Estimation uses CAFE5's Nelder–Mead optimizer on log
rates, multiple deterministic starts, and explicit zero-rate boundary fits.
`--iterations` controls the maximum iterations per start. A finite, converged
best start is required. The optimization table includes starts that failed to
find a finite likelihood. Numerical multistart convergence is not proof of a
global optimum, and no asymptotic confidence interval is implied.

`--likelihood-only` skips family and ancestral output for rate profiling or
sensitivity fits; it retains the results and optimization tables.
Output prefix parent directories must exist. Use the explicit `--gamma-cats`, `--alpha`, `--epsilon`, `--estimate-epsilon`
and `--error-model` options documented in the extension guide. Legacy aliases
`-k`, `-e` and `--zero_root` are not accepted in this mode.

Outputs:

| File | Contents |
|---|---|
| `_results.tsv` | Rates, conditional likelihood, root law, convergence and truncation diagnostics |
| `_optimization.tsv` | Individual starts, zero-rate boundaries, scores and convergence |
| `_families.tsv` | Conditional marginal log likelihood per original family ID |
| `_ancestral.tsv` | Nodewise marginal MAP count, posterior mean, probability of zero, and cap mass |
| `_simulated.tsv` | Simulated observed-family CAFE table |
| `_experimental_pvalues.tsv` | Optional fixed-parameter Monte Carlo tail probabilities |

Internal node IDs (`Node0`, `Node1`, ...) use prefix traversal and are specific
to this input tree. Parent IDs and tip names are included. These are not OrthoFinder
node numbers. Nodewise marginal MAP values need not constitute the joint MAP
ancestral assignment; differences between them are not branch significance tests.

Exit status is 0 for successful checked runs, 1 for invalid input or computational
failure, and 2 for a nonconverged best fit or failed final truncation check.
Inspect the diagnostics even when the process returns 0.

## Simulation and significance

```bash
./build/cafe5 --innovation -t tree.nwk --root-mean 1 \
  --lambda 0.002 --nu 0.001 --simulate 5000 --seed 42 -o synthetic
```

Simulation draws unbounded branch counts directly: a binomial number of surviving
initial lineages with geometric descendant counts, plus gamma–Poisson immigration.
It does not truncate or renormalize a transition row. Root zero can generate
positive descendants. Observed-family simulation rejects all-zero tip vectors
under the same inclusion rule used by the likelihood. Root sampling uses the
specified finite root representation. Seeds reproduce runs within the same
C++ standard-library environment; cross-platform bitwise reproducibility is not
promised.

Optional `--bootstrap B` compares each observed conditional log likelihood with
B independent families generated at the supplied/fitted parameter values. It
reports `(1 + number_of_simulated_scores <= observed_score)/(B+1)`.
Simulations outside the state cap cause an explicit error, never silent rejection.
These probabilities are **experimental fixed-parameter model checks**, not
validated CAFE branch p-values. Fitting parameters to the same families can alter
calibration; this option does not refit parameters in each replicate. Multiple
testing correction and a broader assessment of composite-null calibration remain
necessary for inferential use. This legacy family-level option emits no branch p-values; use the separate
`branch_bootstrap.py` driver for branch tests and their calibration workflow.

## Independent validation

```bash
python3 scripts/innovation/validate.py build/cafe5 validation/numerical
python3 scripts/innovation/recovery.py build/cafe5 validation/recovery
```

The first script compares transition probabilities with an independently built
SciPy generator-matrix exponential, tests boundaries, moments, the semigroup law,
small-tree likelihoods and root/internal posteriors, seed reproducibility,
simulation moments, and rejection of invalid options. The second performs 12
independent parameter-recovery experiments and a fixed-parameter null calibration
experiment. They require NumPy and SciPy only for validation, not for CAFE5 itself.
The GitHub workflow also runs the existing CAFE tests. With recent GCC releases,
the unchanged upstream tests may require
`-DCMAKE_CXX_FLAGS="-Wno-error=deprecated-declarations -Wno-error=stringop-overflow"`
because they promote warnings in upstream code to errors.

An independent R transition check is provided in `scripts/innovation/check_kernel.R`.
Numerical equality is evidence about the implementation, not about biological
model adequacy. See `innovation_validation.md` for the actual validation results,
input provenance, empirical sensitivities and remaining limitations.

## Selecting a compatible HOG table

For the verified source OrthoFinder tree used in the wasp example:

- `N10` spans the MRCA of Vespidae and the two external outgroups; subset its columns
  to the 17 tips of `Vespidae_with_outgroups_dated_primary`.
- `N11` corresponds to Vespidae including Ancistrocerus and matches the 15-tip
  `Vespidae_dated_primary` tree.
- `N13` corresponds to the social-wasp crown only, excluding Ancistrocerus.
- `N12` is the other daughter of N10 and does not contain Vespidae.

Empty columns outside a node-specific HOG's scope must not be treated as measured
absences. Node numbers are run-specific; verify these mappings again after an
OrthoFinder rerun or topology change. The conversion script records source hashes,
counts genes in comma-separated membership cells, preserves family IDs, and removes
only rows all-zero in the selected taxa. It does not substitute RNA-seq counts.

```bash
python3 scripts/innovation/prepare_wasps.py /path/to/evodevo_waspcastes validation/wasps
```

The focused example uses **Vespidae_dated_primary + N11**, as selected by the user.
The named Vespidae-only dated tree was dated separately and has some different
branch lengths from simply pruning the full dated tree; its exact file is used.
