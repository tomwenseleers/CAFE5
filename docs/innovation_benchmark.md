# How well does the extended model fit?

On the example dataset, the birth–death–innovation model with gamma rate
variation and count error fits substantially better than models restricted to
equal duplication and loss rates with no innovation. The improvement remains
large after AIC penalizes the additional fitted parameters.

This page explains the comparison, how the number of gamma categories was
chosen, and how we checked that numerical approximations did not determine the
result. The [example data](../examples/innovation_benchmark/) contain 13,742
observed gene families across 15 species. The largest observed count is 145
copies in one species. **All 13,742 supplied families were fitted together.**
Input sources and preparation are recorded in the example's
[provenance file](../examples/innovation_benchmark/provenance.json).

## Why extend the birth–death model?

The original motivation was the zero-root problem. Duplication can increase the
size of an existing family, but cannot create a copy of a family that is absent.
Under a pure duplication–loss process, a family with zero copies at the root
must therefore remain absent throughout the tree. Adding a copy-independent
innovation rate allows gains from zero. With a root distribution that includes
zero, families can originate after the root and still be analysed on the full
tree.

The extension also allows duplication and loss to have different rates, families
to evolve at different speeds, and observed counts to contain error. These
features address different aspects of the data. The benchmark below evaluates
their combined contribution; it does not attribute the entire improvement to
innovation alone.

For this particular comparison, **every model assumes one copy at the root**.
This keeps the root assumption identical across models. It measures the benefit
of the broader process and error model, including gains after extinction, but
does not test reconstruction of a zero root. To allow root absence in an analysis,
use `--root poisson`, as described in the [usage guide](innovation.md).

## Comparing the models fairly

We fitted three models to exactly the same count table and tree:

1. **Equal-rate birth–death:** duplication and loss share one rate, innovation
   is absent, and all families share the same evolutionary rate. One count-error
   parameter is also estimated, giving two fitted parameters in total.
2. **Equal-rate birth–death with gamma variation:** the same model, with a gamma
   distribution allowing some families to evolve faster than others. Estimating
   its shape adds one parameter, giving three in total.
3. **Birth–death–innovation with gamma variation and count error:** duplication,
   loss and innovation have separate rates. The model also estimates the gamma
   shape and separate error parameters for zero and positive true counts,
   giving six fitted parameters in total.

The first two models use the core evolutionary assumptions of CAFE5's base and
gamma models. We fitted them within this extension so that all three use the
same root assumption and the same likelihood calculation, including conditioning
on families being observed in at least one species. They are not runs of the
unmodified CAFE5 executable: its treatment of root states and likelihood
normalization differs. The [mathematical comparison](innovation_mathematics.md)
explains those differences.

AIC combines goodness of fit with a penalty for estimating more parameters:

$$
\mathrm{AIC}=2\,\mathrm{NLL}+2p,
$$

where NLL is the negative log likelihood and *p* is the number of fitted
parameters. Lower AIC is better. ΔAIC is the difference from the best model in
the table.

| Model | Fitted parameters | NLL | AIC | ΔAIC |
|---|---:|---:|---:|---:|
| Equal-rate birth–death | 2 | 117707.021 | 235418.041 | 50895.080 |
| Equal-rate birth–death with gamma variation | 3 | 95606.972 | 191219.945 | 6696.983 |
| Birth–death–innovation with gamma variation and count error | 6 | 92255.481 | 184522.961 | 0.000 |

Allowing gamma rate variation greatly improves the fit of the equal-rate model.
The full extension improves AIC by a further **6,696.983**, even after the penalty
for its three additional parameters. This is strong relative support for the
extended model among these candidates on this dataset. It does not establish
that the same model will be preferred for every dataset, or identify which
individual addition explains most of the gain.

## Choosing the number of gamma categories

The gamma distribution describes continuous variation in evolutionary rates
among families. Computation approximates that distribution using *K* rate
categories. These categories are determined by a single fitted shape parameter;
they are not *K* independently estimated rates.

The following table shows the verified category-refinement fits for the
six-parameter model:

| Gamma categories K | AIC | ΔAIC |
|---:|---:|---:|
| 11 | 184531.711 | 8.750 |
| 12 | 184525.314 | 2.353 |
| 13 | 184522.961 | 0.000 |
| 14 | 184523.439 | 0.478 |
| 16 | 184527.050 | 4.089 |

We selected **13 categories**, which gave the lowest AIC among these verified
fits. Fourteen categories gave almost the same result. This supports using 13
as a practical numerical approximation here, rather than interpreting the data
as evidence for exactly 13 biological classes of families. All rows estimate
six parameters, so their AIC differences reflect differences in fit. The best
category count should be checked again for a new dataset.

The selected fit estimates duplication λ = 0.01893405, loss μ = 0.03325355,
innovation ν = 0.00026007, gamma shape α = 0.17720017, positive-count error
ε = 0.01045103, and zero-count error ε₀ = 0.05291919. The evolutionary rates
are expressed per unit of time in the supplied tree; the gamma shape and error
parameters are dimensionless. The fitted loss rate is approximately 1.76 times
the duplication rate.

## Keeping large families in the analysis

Large families and large changes in copy number can make likelihood calculations
numerically difficult. The workflow keeps these families in the analysis and
increases the range of possible ancestral counts used in the calculation when
needed. It then checks whether increasing that range further changes the
likelihood appreciably. Small and large families contribute to the same fitted
model and bootstrap analysis.

The upper numerical count limit is called the **cap**. It is a computational
setting, not a threshold for removing observed families. In this benchmark,
the equal-rate gamma model required a cap of 360; the other two models used 180.
Each accepted fit had to change its total NLL by less than 0.01 when the cap was
doubled. The equal-rate gamma fit changed by approximately 0.000290, and the
selected model changed by zero at the reported precision. These numerical
differences are much smaller than the differences between the fitted models.

The workflow also uses multiple optimization starts and independent checks of
likelihood normalization. The selected likelihood was reproduced with the
release binary. These checks provide evidence that the reported improvement is
not an artefact of the count limit, although multiple starts cannot prove that
a global optimum has been found. Available memory, runtime and the workflow's
maximum cap still limit the sizes of problems that can be analysed.

## What the benchmark establishes

The results support the combined extension over the two simpler models under
shared assumptions on this dataset. Establishing how much each feature
contributes would require further comparisons that remove one feature at a
time. Testing performance on independent datasets and simulated data addresses
a different question from comparing AIC on this example.

A lower AIC also does not by itself establish the accuracy of family or branch
p values. Those are calculated using the separate
[refitted-bootstrap procedure](focal_bootstrap_methods.md), and depend on the
fitted model adequately representing the data. They are nominal p values;
no FDR adjustment is applied.

## Reproduce the comparison

Run this command from the repository root after building the executable:

```bash
python3 scripts/innovation/benchmark_models.py build/cafe5 \
  examples/innovation_benchmark/counts.tsv \
  examples/innovation_benchmark/tree.nwk -o benchmark --threads 4
```

This refits the three models in the first table. It increases numerical caps
when checks fail and saves the commands, diagnostics and `comparison.tsv`.
The full comparison can take substantial time. To explore gamma category counts
and then bootstrap a selected model, follow the [usage guide](innovation.md).

The saved [release evidence](validation/release/) includes the
[model comparison](validation/release/comparison.tsv),
[gamma category comparison](validation/release/gamma_resolution.tsv), and
independent numerical checks. It also includes a fixed-parameter calculation
that reproduces the selected likelihood. That calculation verifies a previous
fit; it does not re-estimate parameters or change the original model's
six-parameter AIC penalty.

To verify the AIC arithmetic in R:

```r
z <- read.delim("benchmark/comparison.tsv")
stopifnot(all(abs(z$AIC - (2*z$NLL + 2*z$parameters)) < 1e-7))
z$delta_AIC <- z$AIC - min(z$AIC)
z[order(z$AIC), ]
```
