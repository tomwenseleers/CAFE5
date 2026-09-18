# Refitted bootstrap: protocol and validation

The supported driver is `scripts/innovation/refitted_bootstrap.py`, normally
invoked by `analyze.py`. It supports unequal λ/μ, unscaled innovation, discrete
gamma variation, separate zero-error estimation, and parametric root laws. It
replays the original fit's fixed/estimated parameter choices explicitly.
Nonparametric root files and error tables are not currently accepted by this
bootstrap driver. [Equations and comparison with CAFE5](innovation_mathematics.md)
define the tests and their limitations.

Every simulated dataset has the same number of observed families as the input;
every dataset is refitted. No large family or failed replicate is silently
removed. Caps grow to accommodate simulated counts and must pass a doubled-cap
likelihood check. Failed optimization triggers a documented retry, followed by an
explicit error if unsuccessful. The default performs boundary fits; the advanced
interior-refit option audits initial datasets and falls back near selected
boundaries, and should be disclosed if used.

`branch_p` is the doubled-minimum tail of signed posterior mean change;
`transition_branch_p` is a sensitivity statistic. Output is nominal, without FDR
or familywise correction. Tests cover all branches by default; `--nodes` restricts
the output/reference calculations to specified nodes. A selected branch is
unusual under the fitted global process, not a test of zero events on that branch.

The bootstrap stores independent-dataset Monte Carlo intervals. Their width is
not biological effect uncertainty. Simulation count B times family count n must
not be described as Bn independent model refits. The per-dataset plus-one tail
correction implies finite resolution that increasing B alone cannot overcome.

## Reproducible checks

The GitHub workflow `.github/workflows/innovation.yml` builds the engine and runs
legacy unit tests plus independent transition/likelihood, gamma/error, root-law,
conditional simulation, rare-inclusion and broad-tail normalization, and bootstrap checks. Important entry points:

```bash
./build/tests
python3 scripts/innovation/validate.py build/cafe5 validation/numerical
python3 scripts/innovation/validate_improvements.py build/cafe5 validation/extended
python3 scripts/innovation/validate_rare_likelihood.py build/cafe5 validation/rare
python3 scripts/innovation/validate_workflow.py build/cafe5 validation/workflow
```

The local release checks include 206 legacy test cases (582 assertions),
independent extended-model checks, and an end-to-end two-dataset bootstrap smoke
check including resume and incompatible-input rejection. A two-dataset smoke
check verifies execution and output contracts, **not statistical calibration**.
Existing small simulation/recovery tests likewise cannot establish uniform tail
calibration for every fitted model or data-generating regime. Substantive analyses
should inspect predictive adequacy, convergence, cap stability and Monte Carlo
uncertainty; a better in-sample AIC does not replace these checks.
