# Manuscript comparison: full-data bootstrap protocol

The primary comparison uses the original **N0 HOG definitions** and the dated
**17-species Vespidae-plus-outgroups tree**, as requested on 18 September 2026.
The raw N0 membership source has 24,502 rows. Of these, 12,577 have at least one
member in these 17 species. All 12,577 enter a single joint likelihood, including
1,550 single-species families and 56 families with a maximum-minus-minimum count
above 20. The maximum observed count is 141. No family is removed for inferred
root state or count differential. All-zero rows contain no observation in this
sample and are outside the observed-family sampling distribution.

The model retains equal per-copy duplication and loss, global innovation, discrete
gamma variation on duplication/loss only, and the symmetric one-copy observation
error model. The current comparison conditions on at least one observed copy and
jointly estimates the mean of a zero-inclusive Poisson root law. A profile check
showed that the prototype's fixed mean of 1 was inadequate for this dataset. The
root-law parameter is distinct from each family's posterior ancestral count. Three- and six-category gamma fits are being compared; the selected
category count is held fixed during each model's bootstrap.

## Statistics and thresholds

For observed family y and fitted parameters theta, the family statistic is

    S(y; theta) = log Pr_theta(Y=y | Y is not all zero).

Its p-value is the inclusive lower tail under the fitted population of observed
families. The original branch sensitivity statistic is

    T_b(y; theta) = E_theta[N_child - N_parent | Y=y].

Root states, gamma categories, and annotation error are marginalized in this
expectation. For this sensitivity statistic, the branch p-value is twice the smaller inclusive
tail, capped at 1.

The primary branch comparison uses a transition-surprise statistic closer to the
manuscript's transition-based question. Let i and j be the marginal-posterior MAP
counts at the parent and child, and w_k(y) the family posterior gamma-category
weights. Define

    L_b = sum_k w_k(y) Pr_k(N_child <= j | N_parent=i),
    U_b = sum_k w_k(y) Pr_k(N_child >= j | N_parent=i),
    R_b = -log min(1, 2 min(L_b, U_b)).

The primary branch p-value is the inclusive **upper simulated tail of R_b**.
The quantity inside the logarithm is NOT itself claimed to be a calibrated
p-value: counts and category weights were inferred from the same observations.
That entire inference procedure is repeated on the bootstrap data. The upper
transition tail includes mass beyond the finite numerical cap. Direction and
integer count change use j-i. The requested primary calls require family p <
0.05, calibrated transition-surprise branch p < 0.01, and j-i != 0. Both this
statistic and the posterior-mean sensitivity are reported, rather than choosing
between them based on which reproduces more manuscript results.

These are population-model tail tests. The branch null is not literal absence of
all gains and losses on that branch. In particular, the branch p-value is not the
legacy Viterbi probability. Equal numerical thresholds do not make those two
statistics equivalent. No multiplicity adjustment is used to reproduce the
requested decision rule; BY and Holm values across the two focal branches are
provided as supplementary columns, subject to the bootstrap approximation.

## Pooled, refitted parametric bootstrap

Each replicate simulates an entire dataset with the same number of observed
families at the fitted parameters. All five originally free parameters (lambda,
nu, alpha, epsilon, Poisson root mean) are then refitted with the same root
law family, error model, gamma category count, and ascertainment correction. Repeated
count patterns are likelihood-compressed without changing their multiplicities.
The first optimizer start uses the generating fit. The empirical fit searches
exact-zero lambda, nu, and epsilon faces. To reduce computation, the first two
bootstrap datasets audit those boundary alternatives; subsequent datasets refit
all five parameters in the positive interior. A tenfold drop in lambda, nu, or
epsilon relative to its generating value triggers a full boundary-search fallback.
This is an optimization shortcut, not a Gaussian approximation to the test-statistic
distribution. Boundary audit results are retained. A small-tree comparison of full
and interior fits agreed in NLL to 4e-9. Failed fits are retained and stop the run, rather than being
silently omitted or replaced by a more convenient random seed.

For replicate r with n simulated families, define the lower-tail estimate

    F_r(s) = (1 + sum_i 1[S(Y_ri; theta_hat_r) <= s]) / (n+1).

The reported family p-value is the average of F_r evaluated at the observed
statistic. Transition-surprise upper tails are calculated analogously. For the mean-change
sensitivity, both tails are averaged before taking twice the smaller mean. Because families are exchangeable under this particular global
model, pooling estimates the same marginal bootstrap CDF as taking a single
family index per simulated dataset. It uses more of each expensive refit.
However, the n fitted statistics within a dataset are dependent through their
shared parameter estimates. **B datasets are B independent clusters, not B*n
independent parameter refits.** Monte Carlo standard errors are therefore the
sample standard deviation of the B dataset-level tail estimates divided by
sqrt(B). Approximate pointwise 95% Student-t intervals are reported; they measure
Monte Carlo uncertainty, not biological effect uncertainty or model adequacy.
The +1 smoothing avoids zero estimates but is not claimed to yield an exact
finite-sample p-value for the pooled, fitted-parameter procedure.

This approach reduces the number of complete refits needed for useful tail
resolution without assuming a Gaussian or chi-square distribution for the branch
statistic. The parametric bootstrap itself remains an approximation when nuisance
parameters are estimated. Threshold crossings of the Monte Carlo intervals are
flagged rather than concealed.

## Numerical support and family retention

Simulation uses the unbounded process. The latent numerical count cap is set above
the largest observed count in each simulated dataset, and every fitted dataset is
checked by doubling the cap. A failed numerical tail check causes a larger-cap
refit of that same dataset. It never causes removal of a large family. A practical
resource ceiling can stop a run for investigation; it is not a biological size
filter. The observed data, parameters, random seeds, logs, and fit diagnostics are
retained for auditing.

TE-associated families remain in the fit, matching the manuscript's fitting
policy. Both all-family calls and TE-excluded biological summaries are supplied.
Direct member annotations and the existing same-N0 annotations are used in the
comparison. The earlier N11 pilot is a sensitivity analysis only: 22 manuscript
focal events overlap 136 N11 families and are not one-to-one family identities.

## Why pooling preserves the target CDF

Conditional on the empirical generating estimate theta_hat, the simulated
families in a replicate are exchangeable. The estimator theta_hat_r is a
permutation-invariant function of that entire replicate. Consequently the fitted
statistics S(Y_ri; theta_hat_r) are exchangeable even though they are dependent.
For any fixed threshold s,

    E[(1/n) sum_i 1{S(Y_ri; theta_hat_r) <= s}]
      = Pr{S(Y_r1; theta_hat_r) <= s}.

Thus pooling is not a different null hypothesis. Its within-dataset variance is

    Var(F_r) = Var(I_r1)/n + ((n-1)/n) Cov(I_r1, I_r2),

before the small +1 smoothing transformation. Treating B*n observations as
independent drops the covariance term and can understate Monte Carlo uncertainty.
The implemented between-dataset variance estimate retains it. This argument does
not make plug-in nuisance estimation exact, remove model misspecification, or
justify treating a small simulation as proof of calibration at every parameter.

## Input identity check

The 17-tip Newick text matches the archived manuscript CAFE input exactly. Counts
for all 11,025 originally prepared families also match exactly. The new input adds
1,550 single-species families and two multispecies families (`N0.HOG0000000` and
`N0.HOG0000119`) absent from that preparation. It retains the old large-family
track within the joint fit. In the innovation engine's prefix-order numbering,
`Node3` corresponds to manuscript node 20 (social Vespidae) and `Node13` to node
25 (Vespinae). These identities are verified from descendant species sets.

## Interpreting the two branch statistics

The two p-values answer different questions. The posterior-mean test asks whether
the inferred signed number of copies gained or lost is extreme among exchangeable
families from the global model. It need not be extreme relative to that family's
own inferred turnover category and ancestral size. The transition-surprise test
uses the inferred starting size and category weights, but relies on marginal-MAP
ancestral states. Marginal modes at two nodes need not be the jointly most likely
pair of states. A rare transition between those modes can therefore coexist with
substantial posterior uncertainty about the direction of change.

For that reason, the manuscript comparison reports **both tests**, with explicitly
labelled tables, and independently calculated joint parent/child posterior
diagnostics for the published and newly flagged HOGs. Posterior expansion
probabilities are not p-values. Neither the choice of statistic nor the root law
is justified by whether it recovers more of the manuscript's calls. A mismatch in
predictive occupancy distributions limits biological interpretation even when
simulation-based calibration under the assumed model passes.

TE annotation uses the complete archived 23-HOG exclusion list, the archived
root-zero HOG annotations, and direct member evidence; it does not rely only on
the much shorter list of TE HOGs significant at the two focal nodes. Newly flagged
families are screened for explicit TE terms as well. Every such flag is retained
in the all-family output, and no TE flag removes a family from the likelihood.
