# Birth–death–innovation inference and its relationship to CAFE5

This document specifies the implemented likelihood and tests. The reference for
standard CAFE5 is
[Mendes et al. (2020), *CAFE 5 models variation in evolutionary rates among gene
families*](https://doi.org/10.1093/bioinformatics/btaa1022), with implementation
in `src/base_model.cpp`, `src/gamma_core.cpp`, `src/probability.cpp`,
`src/execute.cpp` and `src/gene_family_reconstructor.cpp`. The extension is in
`src/innovation.cpp`; the supported inferential driver is
`scripts/innovation/refitted_bootstrap.py`.

## 1. Count process

Let $`N(t)`$ be the true copy count of one family along an edge of length $`t`$.
Conditional on the family rate multiplier $`r`$, the continuous-time Markov generator is

```math
q_{n,n+1}=r\lambda n+\nu,\qquad
q_{n,n-1}=r\mu n\;(n\geq1),\qquad
q_{n,n}=-(r\lambda n+r\mu n+\nu).
```

All other instantaneous rates vanish. The transition probability is
$`P_{ab}(t;r,\theta)=\Pr(N(t)=b\mid N(0)=a,r,\theta)`$.
The implementation evaluates analytic birth–death–immigration transition
probabilities, including the equal-rate limit, rather than estimating a separate
parameter for every matrix entry. Counts at internal nodes are latent.

For standard CAFE5, $`\nu=0`$ and $`\lambda=\mu`$; hence state zero is absorbing:
$`P_{0b}(t)=0`$ for $`b>0`$. With $`\nu>0`$, $`q_{0,1}=\nu`$, making acquisitions after
absence possible. This resolves the process-level zero-state restriction; the
chosen root distribution must also assign positive mass to zero if root absence
is part of the scientific model.

For $`d=r(\lambda-\mu)`$, the conditional expectation satisfies

```math
\frac{d}{dt}E[N(t)]=dE[N(t)]+\nu,
\qquad
E[N(t)\mid N(0)=a]=
\begin{cases}
ae^{dt}+\nu(e^{dt}-1)/d,&d\ne0,\\[2pt]
a+\nu t,&d=0.
\end{cases}
```

Thus unequal rates allow a systematic tendency toward net loss or growth,
whereas the critical duplication/loss process has no such drift. Innovation is
constant per family per time unit, whereas duplication is proportional to count.

## 2. Gamma heterogeneity and count error

The continuous family multiplier has mean-one distribution
$`r\sim\mathrm{Gamma}(\alpha,\text{rate}=\alpha)`$. Numerical integration uses
K equally weighted categories: $`w_k=1/K`$ and $`r_k`$ the midpoint quantile of the kth
probability interval, rescaled so the discrete mean is one, as in the active
median-discretization path of CAFE's gamma routine. A family retains its
category throughout the tree. Only duplication and loss are scaled; innovation
remains $`\nu`$. For K=1, $`r_1=1`$ and $`\alpha`$ is not estimated.

For a true positive tip count $`n`$, the observation law is

```math
E(y\mid n)=
\begin{cases}
\varepsilon,&y=n-1\text{ or }n+1,\\[2pt]
1-2\varepsilon,&y=n,\\[2pt]
0,&\text{otherwise},
\end{cases}
\qquad n\geq1.
```

At zero, $`E(1\mid0)=\varepsilon_0`$ and
$`E(0\mid0)=1-\varepsilon_0`$. This distinguishes zero-to-one annotation error from
one-step error around positive counts. It is a restricted observation model,
not an unrestricted annotation correction. The optimizer constrains $`0\leq\varepsilon<1/2`$,
$`0\leq\varepsilon_0<1`$ and (when fitted) $`0.05\leq\alpha\leq100`$.
An alpha estimate at a bound requires scrutiny. The same law applies to every tip;
lineage-specific annotation or ascertainment effects are not modeled.

## 3. Root law, pruning and ascertainment

The simple preset has $`\pi(a)=\mathbf{1}\{a=1\}`$. The Poisson option uses
$`\pi(a)=e^{-\rho}\rho^a/a!`$, estimating $`\rho`$ and including root zero.
The low-level interface also offers hurdle laws:
$`\pi(0)=z`$ and $`N_{root}\mid N_{root}>0=1+X`$, with $`X`$ Poisson or negative
binomial. The mean parameter describes $`X`$, not the unconditional root mean.

For family i with observed vector $`y_i`$, the pruning recursion is

```math
F_{v,k}(a)=
\begin{cases}
E(y_{iv}\mid a),&v\text{ a tip},\\[2pt]
\prod_{c\in children(v)}\sum_bP_{ab}(t_c;r_k,\theta)F_{c,k}(b),&\text{otherwise}.
\end{cases}
```

The unconditional likelihood is a **sum** over root states and categories:

```math
L_i(\theta)=\sum_k w_k\sum_a\pi_\theta(a)F_{root,k}(a).
```

The likelihood treats families as independent draws from the common model.
The data contain only observed families, $`A=\{\sum_vY_v>0\}`$. Therefore

```math
L_i^{+}(\theta)=\frac{L_i(\theta)}{\Pr_\theta(A)},\qquad
\ell(\theta)=\sum_{i=1}^n\log L_i^{+}(\theta).
```

Mix categories **before dividing** by the mixture inclusion probability.
Conditioning each category and retaining its original weight would generally
be a different likelihood. Observation error is included in $`\Pr(A)`$.

For ordinary inclusion probabilities, the implementation evaluates
$`\log\Pr(A)=\log[-\mathrm{expm1}(\log\Pr(Y=0))]`$. This avoids
enumerating an arbitrarily broad tail of possible positive tip counts. When
$`\log\Pr(Y=0)\geq-0.01`$, subtraction becomes progressively ill-conditioned,
so the implementation instead propagates both zero-observation mass Z and
positive-observation mass U. Combining
an accumulated subtree with a child whose transition-integrated masses are
$`(z,u)`$ uses $`U_{new}=U(z+u)+Zu`$ and $`Z_{new}=Zz`$. Root and category sums of U
give $`\Pr(A)`$ directly. Both formulas represent the same infinite-state
probability; finite-cap error is controlled by the doubled-cap check. Analytic
tests cover both sides of the numerical switch, probabilities down to
approximately $`10^{-120}`$, and a broad critical-process tail with count 145.

Internal counts use a finite numerical cap; transition rows are not renormalized
after truncation. A fit must pass a doubled-cap likelihood check
$`|\mathrm{NLL}_{M}-\mathrm{NLL}_{2M}|<0.01`$ and root omitted-mass checks. The
workflow enlarges the cap when necessary and retains every supplied family.
This is an explicit numerical tolerance, not a proof of zero truncation error.

### Contrast with the legacy fitted objective

In this checkout, the base CAFE objective profiles over the positive root state,
using the largest root-prior-weighted pruning likelihood. The gamma objective
profiles within each category and then sums across categories. Schematically,

```math
\widetilde{L}_{i,base}=\max_{a\geq1}\{\pi(a)F_i(a)\},\qquad
\widetilde{L}_{i,gamma}=\sum_k w_k\max_{a\geq1}\{\pi(a)F_{i,k}(a)\}.
```

These are not the root-integrated, observed-family-conditioned probabilities
above. Consequently raw legacy objective values and extension AIC values should
not be placed in one AIC ranking. The benchmark instead refits CAFE-like process
restrictions within a common normalized likelihood, root law and data universe.

## 4. Family p values

### Standard CAFE calculation

`compute_pvalues` simulates 1,000 reference families for each candidate positive
root count under the fitted base duplication/loss rate. For each root it compares
the observed conditional pruning likelihood with simulated likelihoods and takes
the largest tail probability over the root counts searched. Schematically,

```math
p_i^{legacy}=\max_{a\in\mathcal{R}_{i}}
\widehat{\Pr}_a\{F(Y;a)\leq F(y_i;a)\}.
```

The code searches a family-size-dependent root range; its finite-sample index and
tie rules are defined in `pvalue`/`find_best_pvalue`. These reference simulations
do not refit the dataset. In `execute.cpp`, the gamma analysis also calls this
routine with the original/base lambda, not the full fitted gamma mixture.
The p-value pruning call does not pass the fitted error model, so it should not
be described as a full gamma-plus-error refitted calibration.

### New refitted pooled reference distribution

Let $`\hat{\theta}`$ be the selected model's fit. Generate B complete independent
observed-family datasets, each with n families, from $`\hat{\theta}`$. Refit the
**same free parameters** on each dataset to obtain $`\hat{\theta}_b^{*}`$. Compute

```math
T_i=\log L_i^{+}(\hat{\theta}),\qquad
T_{bj}^{*}=\log L^{+}(Y_{bj}^{*};\hat{\theta}_b^{*}).
```

Within each replicate form the corrected empirical lower tail, then average:

```math
f_{bi}=\frac{1+\sum_{j=1}^n\mathbf{1}\{T_{bj}^{*}\leq T_i\}}{n+1},
\qquad \hat{p}_i=\frac{1}{B}\sum_{b=1}^B f_{bi}.
```

Pooling assumes families are exchangeable under the fitted global model; this is
a model-based approximation, not an exact test conditional on each family's root
size or gamma category. The plus-one correction prevents zero estimates; it does
not make the refitted plug-in bootstrap an exact finite-sample test. Model
selection is held fixed, so its uncertainty is not included. HOG reconstruction
can impose additional ascertainment beyond presence in at least one sampled
species; only this latter observation condition is modeled here.

## 5. Branch p values and direction

### Standard CAFE calculation

For reconstructed parent/child states a,b and transition matrix P, the legacy
branch calculation is a probability-ordering mid-p:

```math
p_{\mathrm{edge}}^{\mathrm{legacy}}
=\sum_{m\,:\,P_{am}\lt P_{ab}} P_{am}
+\frac{1}{2}\sum_{m\,:\,P_{am}=P_{ab}} P_{am}.
```

The implemented sum has a finite count range. It conditions on reconstructed
ancestral counts rather than integrating their uncertainty or refitting simulated
datasets. `execute.cpp` calculates these branch values only for families passing
the family threshold. Neither stage automatically applies FDR correction.

### New primary branch statistic

For each branch e, use the signed posterior mean change

```math
D_{ie}=E_{\hat{\theta}}[N_{child(e)}-N_{parent(e)}\mid y_i].
```

The posterior averages over parent/child states and family gamma categories,
using all tip data and the observation model. For each refitted simulated family,
compute $`D_{bje}^{*}`$ by exactly the same procedure. At each branch,

```math
u_{bie}=\frac{1+\sum_j\mathbf{1}\{D_{bje}^{*}\geq D_{ie}\}}{n+1},\qquad
l_{bie}=\frac{1+\sum_j\mathbf{1}\{D_{bje}^{*}\leq D_{ie}\}}{n+1},
```

```math
\hat{p}_{ie}=\min\{1,\;2\min(\bar{u}_{ie},\bar{l}_{ie})\}.
```

Expansion/contraction follows the sign of $`D_{ie}`$, not a separately rounded
ancestral MAP change. The null is **the fitted global evolution model for an
exchangeable observed family**, not “no change on this edge”. A significant
contraction means unusually negative change relative to that null, not simply
posterior evidence that at least one loss occurred.

For the transition sensitivity, let $`a^{*},b^{*}`$ be marginal posterior modes and
$`w_{ik}=\Pr(k\mid y_i)`$. Compute

```math
H^{-}_{ie}=\sum_k w_{ik}\sum_{m\leq b^{*}}P_{a^{*}m}(t_e;r_k),\qquad
H^{+}_{ie}=\sum_k w_{ik}\sum_{m\geq b^{*}}P_{a^{*}m}(t_e;r_k),
```

and $`S_{ie}=-\log\min\{1,2\min(H^{-}_{ie},H^{+}_{ie})\}`$, with a machine-minimum
floor inside the logarithm. The reported sensitivity p value is the pooled
refitted bootstrap upper tail of S, using the same plus-one convention. This
statistic conditions on marginal modes and is distinct from posterior mean
change; do not mix its calls with the primary directional report.

## 6. Monte Carlo uncertainty, selection and limitations

The independent Monte Carlo units are B simulated **datasets**, not Bn independent
refits. For a one-sided tail estimate $`\bar{x}`$, report the approximate interval

```math
\bar{x}\pm t_{B-1,0.975}\frac{s_x}{\sqrt{B}},
```

clipped to [0,1]. For doubled-minimum branch tails the implementation propagates
the two marginal endpoint intervals through the same min/doubling operation.
These are approximate pointwise Monte Carlo intervals, not guaranteed joint
coverage, biological effect intervals, or a substitute for model adequacy checks.
Inspect threshold-crossing intervals before making definitive borderline calls.

The estimator has family floor $`1/(n+1)`$ and two-sided branch floor $`2/(n+1)`$;
increasing B reduces Monte Carlo variability but does not lower these floors.
Nominal reporting uses family p < .05 with branch p < .01, and .05 as sensitivity.
There is **no FDR adjustment**. Applying both thresholds is not a theorem of
experiment-wide error control. Calibration depends on the selected model,
exchangeability, convergence and numerical accuracy. Better AIC alone establishes
neither correct biological annotation nor calibrated false-positive rates under
model misspecification.

For model selection,

```math
\mathrm{AIC}=2p-2\ell(\hat{\theta}).
```

The root-one BDI/gamma/error model has p=6 for K>1 and p=5 for K=1; estimating a
Poisson root mean adds one. K controls deterministic quadrature resolution, not K
free rate parameters. Boundary estimates and model search make ordinary AIC an
approximation; neighboring K with tiny differences should be treated as similar.

## Reproducing nominal selection in R

```r
x <- read.delim("innovation_results/bootstrap/branch_tests.tsv",
                check.names = FALSE)
x$selected <- with(x, family_p < 0.05 & branch_p < 0.01 &
                        posterior_mean_change != 0)
with(subset(x, selected), table(Node, direction))
# Sensitivity analysis, preserving the same family threshold:
with(subset(x, family_p < 0.05 & branch_p < 0.05 &
               posterior_mean_change != 0), table(Node, direction))
```
