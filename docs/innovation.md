# Innovation model: installation and use

This branch adds a birth–death–innovation (BDI) model to CAFE5. It jointly fits
unequal duplication and loss rates, a global innovation rate, gamma variation
among families, and a two-parameter count-error law. A separate, reproducible
whole-dataset bootstrap refits the model and reports **nominal** family and branch
p values. The original CAFE5 interface remains available without `--innovation`.

- [Mathematical model and comparison with CAFE5](innovation_mathematics.md)
- [Comparable-likelihood fit benchmark](innovation_benchmark.md)
- [Bootstrap interpretation and numerical validation](focal_bootstrap_methods.md)

## Install (Linux or Ubuntu under WSL2)

```bash
sudo apt-get install build-essential cmake zlib1g-dev python3-numpy python3-scipy
git clone --branch innovation-bdi https://github.com/tomwenseleers/CAFE5.git
cd CAFE5
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-Wno-error=deprecated-declarations -Wno-error=stringop-overflow"
cmake --build build --target cafe5 tests -j 2
./build/tests
```

## Two-input analysis

```bash
python3 scripts/innovation/analyze.py counts.tsv tree.nwk
```

The table has the usual CAFE columns: description, unique family ID, then integer
counts for each tree tip. Supply a rooted tree with nonnegative branch lengths in
consistent time units. All-zero families cannot be included in an analysis
conditioned on being observed. No other family is dropped because of its size or
inferred ancestral count. Apply any biological exclusions (for example, a
specified transposable-element screen) **before** fitting. The software does not
infer TE annotations or exclude families automatically.

The preset estimates six parameters: duplication λ, loss μ, innovation ν,
gamma shape α, positive-count error ε, and zero-count error ε₀. It uses 13
equal-weight gamma categories, **fixes the root count to one**, and runs 96 full
bootstrap datasets. Gamma multiplies duplication/loss only. Rates are global
across the tree. Root count one is an explicit clade-origin HOG assumption, not
an inference made from arbitrary input. Thirteen categories were useful in the
illustrative benchmark; they are not universally optimal.

For a family universe that permits absence at the root:

```bash
python3 scripts/innovation/analyze.py counts.tsv tree.nwk \
  --root poisson -o poisson_results
```

This estimates a Poisson root mean and allows root zero, adding one fitted
parameter. Innovation can then create copies after ancestral absence. An
innovation process also allows reappearance after extinction on any branch in
the root-one model. Innovation represents acquisition of membership in a family;
it is not a direct molecular mechanism or a per-copy duplication event.

## Resource controls and selected branches

```bash
python3 scripts/innovation/analyze.py counts.tsv tree.nwk \
  -o results --workers 3 --threads 4 --replicates 96
# First inspect fit_branch_statistics.tsv to identify nodes.
python3 scripts/innovation/analyze.py counts.tsv tree.nwk \
  -o fit_only --fit-only
python3 scripts/innovation/analyze.py counts.tsv tree.nwk \
  -o fit_only --nodes Node2 Node12
```

Internal nodes are numbered in the engine traversal; use the `parent` and `Node`
columns in `fit_branch_statistics.tsv` to map the supplied topology. Workers × threads bounds parallel
fit threads. Refits can take hours or days. Interrupting and rerunning the same
command resumes completed datasets. Input/model/binary hashes prevent reusing
incompatible results; use a new output directory after a software/model change.
Increasing `--replicates` reuses existing compatible bootstrap datasets.

## Results

- `analysis.json`: input/model/software provenance.
- `fit_results.tsv`, `fit_complete.json`: parameters, numerical diagnostics, AIC.
- `fit_families.tsv`, `fit_branch_statistics.tsv`: likelihoods and ancestral statistics.
- `bootstrap/family_tests.tsv`: one nominal family p value and Monte Carlo interval per HOG.
- `bootstrap/branch_tests.tsv`: nominal `family_p`, primary `branch_p`, signed
  posterior mean change, direction, Monte Carlo intervals and selection flags.
- `bootstrap/method.json`, `replicate_*/`: simulation/refit provenance and logs.

`selected_raw_thresholds` uses **family p < 0.05 and branch p < 0.01**;
`selected_branch_005` gives the branch p < 0.05 sensitivity. Direction is the sign
of posterior mean child-minus-parent count. `transition_branch_p` is a distinct
sensitivity statistic, not the primary branch p value. No FDR or familywise
correction is applied. Nominal thresholds do not imply experiment-wide error control.

## Gamma resolution and advanced use

```bash
for k in 1 4 8 12 13 14 16; do
  python3 scripts/innovation/analyze.py counts.tsv tree.nwk \
    --gamma-cats "$k" --fit-only -o "fit_k${k}"
done
```

Compare AIC in each `fit_complete.json` and inspect convergence/cap diagnostics.
K > 1 estimates one gamma-shape parameter, **not K independent rates**. Small AIC
differences between neighboring K describe quadrature sensitivity; do not
interpret them as distinct biological rate classes. The script does not perform
automatic model selection. Bootstrap the selected configuration in its original
output directory by omitting `--fit-only`.

The low-level engine remains available for fixed parameters and root-law
sensitivity analyses (`build/cafe5 --innovation --help`). To reproduce the preset
explicitly:

```bash
build/cafe5 --innovation -i counts.tsv -t tree.nwk -o fit \
  --root-family hurdle-poisson --root-mean 0 --root-zero 0 \
  --gamma-cats 13 --estimate-mu --estimate-epsilon --estimate-epsilon-zero \
  --max-count 180 --iterations 3600 --starts 3 --threads 4
python3 scripts/innovation/refitted_bootstrap.py build/cafe5 \
  --tree tree.nwk --observed fit --output bootstrap \
  --replicates 96 --workers 3 --threads 4 --minimum-cap 180
```

Here 180 is an example numerical cap, not a universal setting. The high-level
workflow chooses and enlarges its cap automatically; the low-level command
requires an appropriate user-supplied cap. Failed refits are retained for diagnosis,
never omitted from calibration. Finite memory/runtime and the workflow cap limit
still apply: support for large families is not a promise of unlimited count sizes.

## Capabilities

| Feature | Innovation workflow | Original CAFE mode |
|---|---|---|
| Unequal duplication and loss | Yes | Equal in CAFE5's fitted birth–death model |
| Innovation from zero | Yes | No |
| Gamma family variation | Yes, λ/μ only | Yes |
| Count error | Estimated ε and ε₀; low-level error tables also available | Legacy error model |
| Branch significance | Refitted bootstrap | Legacy reconstruction/transition calculation |
| Separate evolutionary rate per lineage | Not implemented | Legacy rate-group interface |
| Automatic FDR adjustment | No | No |

Branch significance is fully supported by this workflow. **Branch-specific rate
parameters are a separate model extension and remain unsupported in innovation
mode.** The legacy significance formulas are not substituted into the BDI model;
their different null/reference calculations are documented mathematically.

## Original CAFE5 commands

The original interface is selected by omitting `--innovation`:

```bash
build/cafe5 -i counts.tsv -t tree.nwk -o legacy_base
build/cafe5 -i counts.tsv -t tree.nwk -k 13 -o legacy_gamma
```

These invoke the original root treatment and significance procedures, not the
new refitted bootstrap. Their raw objective values are not directly comparable
with the normalized BDI AIC; use the matched-restriction benchmark instead.
