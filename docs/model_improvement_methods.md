# Optional model extensions for the N11 adequacy experiment

These are experimental options. They do not change the standard CAFE5 mode or
the default equal-rate innovation model. The experiment plan is in
[model_improvement_plan.md](model_improvement_plan.md).

## Evolutionary process

`--mu M` fixes the per-copy loss rate independently of duplication lambda;
`--estimate-mu` estimates it. Omission keeps mu=lambda. Gamma multiplies both
per-copy rates by the same family multiplier g, leaving nu global.

For one ancestor, write the birth/death PGF as

    G(z) = p0 + s*a*z/(1-q*z), a=1-q, s=1-p0.

With r=lambda-mu and w=expm1(r*t)/r (w=t at r=0),
q=lambda*w/(1+lambda*w), p0=mu*w/(1+lambda*w), and
s=exp(r*t)/(1+lambda*w). Immigration has a negative-binomial PGF
[a/(1-q*z)]^(nu/lambda). For lambda=0, surviving ancestors are binomial
and immigrants are Poisson with mean nu*(1-exp(-mu*t))/mu. Their convolution
is evaluated with a positive O(K^2) recurrence. Large positive r*t uses an
exponentially rescaled formula. Tiny per-copy rates require log1p expressions:
innovation is not gamma-scaled, so cancellation in log(a)/lambda is consequential.
Finite transition rows are never renormalized. Simulation samples the unbounded
process, independently of the latent numerical cap.

## Observation law

`--epsilon-zero E` specifies P(Y=1 | N=0); omission ties it to epsilon.
`--estimate-epsilon-zero` estimates this probability, with `--epsilon-zero`
optionally supplying its initial value. Positive true counts retain probabilities
(epsilon, 1-2*epsilon, epsilon) for one-copy undercount, correct count and overcount.
These options conflict with a fixed error-model file. Likelihood and simulation
use exactly the same emission probabilities, including ascertainment.

## Root law

`--root-family poisson` is the default. Alternatives are `hurdle-poisson` and
`hurdle-nb`: P(root=0)=root-zero and root | root>0 = 1+Z, with Z Poisson or NB.
For hurdle laws `--root-mean` denotes E[Z], **not the unconditional root mean**.
The unconditional mean is (1-root-zero)*(1+root-mean). For NB, `--root-shape`
is the dispersion/size; Var(Z)=E[Z]+E[Z]^2/root-shape. `--estimate-root-mean`,
`--estimate-root-zero` and `--estimate-root-shape` estimate the corresponding
parameters. Root shape estimation requires hurdle-nb. Root zero is estimated
on the interior and exact-zero sensitivity fits must be performed separately.
The report includes P_root_zero to distinguish actual zero probability from
an unused hurdle parameter in a Poisson fit.

All parametric root laws are subject to an omitted-tail tolerance of 1e-8 and
an independent doubled-cap likelihood check. The NB parameterization is not a
stationary law imposed on the evolutionary process.

## Identifiability and boundaries

The standard exact-zero audit covers estimated lambda, mu, nu and positive-count
epsilon. It does not exhaustively profile the additional hurdle/error parameters;
explicit fixed-zero and profile fits are required before final inference. The
screening runner deliberately skips boundary fits and reports single-start
results as screening only.

If nu=0 and zero-to-one error=0, a root-zero family cannot be observed. Then
P(Y=y)=(1-pi0)*P(Y=y | root>0) for nonzero y, and the factor (1-pi0) cancels
from the observed-family likelihood. Thus pi0 is exactly unidentifiable in that
limit. No optimizer convergence flag can fix this; near the limit it can also
be weakly identified. Profiles and simulation recovery must distinguish
computational convergence from parameter identification.

## Reproducibility and significance

The validation scripts compare kernels and complete small-tree likelihoods with
independent large-state matrix exponentials, including pure-death, critical,
pure-immigration, and extremely small gamma-scaled rate limits. Simulation
frequencies and one-parameter-at-a-time recovery checks are supplied; these do
not establish joint identification on the wasp tree.

The old pooled-bootstrap and independent posterior scripts explicitly reject
these new parameterizations until their refitting logic is generalized and
validated. This prevents silently simulating or refitting an equal-rate model
when the input fit used unequal rates or a different root/error law. The new
predictive-count script does replay every model option. Predictive simulations
are descriptive checks, not calibrated branch significance tests.
