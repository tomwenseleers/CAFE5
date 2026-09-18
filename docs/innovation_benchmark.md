# Evidence from a common-likelihood benchmark

The [frozen example](../examples/innovation_benchmark/) contains 13,742 observed
families across 15 tips, with maximum observed count 145. All supplied families
are fitted together. The snapshot includes a specified preliminary exclusion
screen, but is **not certified TE-free**; this is a statistical/software example,
not a biological result or annotation-enrichment analysis.

## What is being compared

Every model below uses the **same count table, tree, root count one,
observed-family conditioning and normalized likelihood implementation**:

1. Critical base process: λ = μ, ν = 0, K = 1, tied error ε₀ = ε (p = 2).
2. Critical gamma process: λ = μ, ν = 0, K = 13, tied error ε₀ = ε (p = 3).
3. BDI/gamma/error: separate λ and μ, estimated ν, K = 13, separate estimated ε₀
   and ε (p = 6).

The first two impose the core evolutionary restrictions of CAFE5's base/gamma
models. They are **not raw runs of the upstream CAFE objective**, which profiles
root states and uses a different normalization. This controlled comparison
isolates the benefit of the broader model under shared statistical assumptions;
it does not prove that every possible legacy CAFE configuration will have worse
predictive performance. See [the mathematical comparison](innovation_mathematics.md).

| Model | Fitted parameters | NLL | AIC | ΔAIC |
|---|---:|---:|---:|---:|
| Critical base | 2 | 117707.021 | 235418.041 | 50895.080 |
| Critical gamma | 3 | 95606.972 | 191219.945 | 6696.983 |
| BDI gamma error | 6 | 92255.481 | 184522.961 | 0.000 |

The BDI extension improves AIC by 6697.0 relative to the matched critical gamma model despite estimating three additional parameters. The critical gamma fit uses a cap of 360; the other fits use 180. Cap differences are recorded in the [machine-readable table](validation/release/comparison.tsv).


The comparison uses AIC = 2 NLL + 2p, including every fitted rate/error parameter.
AIC assesses relative in-sample support among these candidates. The improvement
belongs to the combined extension; it cannot all be attributed to innovation
without additional matched ablation fits. Optimization uses multiple starts and
numerical checks, but is not a proof of finding a global optimum.

## Gamma resolution

For the six-parameter root-one BDI/gamma/error model:

| Categories K | AIC | ΔAIC |
|---:|---:|---:|
| 11 | 184531.711 | 8.750 |
| 12 | 184525.314 | 2.353 |
| 13 | 184522.961 | 0 |
| 14 | 184523.439 | 0.478 |
| 16 | 184527.050 | 4.089 |

K=13 was selected from these verified candidates. K=14 is almost equally
supported; the result should not be described as strong evidence for exactly
13 biological classes. K changes the discrete-gamma approximation, while p stays
six. This is a selected, documented preset, not a universal optimum across
family definitions, trees or root laws.

The selected rates are λ = 0.01893405, μ = 0.03325355, ν = 0.00026007,
α = 0.17720017, ε = 0.01045103 and ε₀ = 0.05291919, in the supplied tree's time
units. The loss rate is approximately 1.76 times the duplication rate.

## Numerical and scientific rationale

The extension permits innovation when the ancestral count is zero, including
reappearance after extinction. A Poisson-root analysis also permits absence at
the root; the root-one benchmark itself does not demonstrate root-zero inference.
All count sizes remain in a common likelihood and bootstrap. Adaptive latent
caps and explicit convergence checks replace exclusion based merely on size.
Original CAFE can also process large families; this fork does not claim that
large-family analysis is categorically impossible upstream. Its contribution is
an explicit retention and numerical-validation workflow alongside the broader
process and observation model.

Rare-inclusion and broad-tail normalizations have independent analytic checks.
The selected model's likelihood is reproduced by the release binary, and its
cap-doubling difference is zero at reported precision. All accepted benchmark
fits must pass the same absolute NLL-difference tolerance of 0.01.

Good relative fit does not establish absolute predictive adequacy, correct
annotation, or calibrated false-positive rates under model misspecification.
HOG ascertainment, species-specific annotation differences and model-selection
uncertainty remain limitations. The nominal bootstrap tests retain these model
assumptions; no FDR adjustment is applied.

## Reproduce

```bash
python3 scripts/innovation/benchmark_models.py build/cafe5 \
  examples/innovation_benchmark/counts.tsv \
  examples/innovation_benchmark/tree.nwk -o benchmark --threads 4
```

This refits the three candidates, enlarging failed numerical caps and writing
commands, fit diagnostics and `comparison.tsv`. It can take substantial time.
For the selected-model workflow and bootstrap, use [the two-input interface](innovation.md).
The machine-readable [release evidence](validation/release/) includes AIC and
gamma-resolution tables, parameter score replay and independent numerical checks.
The fixed-parameter replay is a score verification, not an additional parameter
fit; its estimated-parameter flags should not be used to count the six parameters
of the original fitted model.

To verify the AIC arithmetic in R:

```r
z <- read.delim("benchmark/comparison.tsv")
stopifnot(all(abs(z$AIC - (2*z$NLL + 2*z$parameters)) < 1e-7))
z$delta_AIC <- z$AIC - min(z$AIC)
z[order(z$AIC), ]
```
