# N11 model-improvement comparison

All 13,836 observed families are retained on the matching 15-tip Vespidae tree. These are model-development and predictive checks, not new calibrated branch significance results. Convergence alone is not evidence of model adequacy. Failed or unfinished candidates remain visible.

| Model | Status | NLL | Loss / duplication | Single-species % predicted | One copy in every species % predicted |
|---|---|---:|---:|---:|---:|
| baseline | converged_archived_baseline | 97857.587 | 1.000 | 31.028 | 15.068 |
| asymmetric | converged | 97290.305 | 1.488 | 29.103 | 18.096 |
| combined_nb_zero_free | converged | 93627.799 | 1.542 | 13.893 | 27.428 |
| combined_poisson_zero0 | converged | 94876.532 | 1.952 | 16.016 | 24.868 |
| gamma12 | converged | 97850.453 | 1.000 | 30.249 | 15.153 |
| gamma24 | converged | 97796.913 | 1.000 | 29.956 | 15.653 |
| hurdle_nb | converged | 95046.092 | 1.000 | 20.194 | 20.972 |
| hurdle_poisson | converged | 95046.092 | 1.000 | 20.194 | 20.972 |
| zero_error0 | converged | 98240.094 | 1.000 | 24.832 | 17.462 |
| zero_error001 | converged | 98101.523 | 1.000 | 25.894 | 17.395 |
| zero_error_free | converged | 97820.930 | 1.000 | 29.864 | 15.821 |
| root01_asymmetric | converged | 94201.925 | 1.596 | 18.199 | 23.551 |
| root01_zero_error0 | converged | 96411.571 | 1.000 | 18.734 | 20.199 |
| root01 | converged | 95046.092 | 1.000 | 19.868 | 21.071 |
| asymmetric_second_start | converged_boundary_caution | 96674.225 | 2.070 | 15.659 | 32.235 |
| combined_second_start | converged_boundary_caution | 93627.799 | 1.542 | 13.785 | 27.431 |

Observed single-species fraction: **12.764%**. Observed exactly-one-copy-in-every-species fraction: **29.604%**. For hurdle roots, root_mean is the mean number of copies above one conditional on a positive root, not the unconditional mean. Root01 fixes that excess to zero and estimates the root-zero probability.

The predictive tables include mean occupancy, mean count, count tails and per-species checks. Simulation ranges are descriptive replicated-dataset ranges, not confidence intervals for evolutionary changes. Extreme-tail disagreement alone need not invalidate a background model: genuine family-specific changes can be outliers. The widespread occupancy and conserved-copy discrepancies are the main adequacy concern.

## Held-out scores

Models are re-fitted to 11,166 training families and scored without parameter re-estimation on 2,670 test families. Original OG identities define the split, keeping related N11 subfamilies together. These validation scores guide model selection; they are not an unbiased final performance estimate after selecting a model. The root01 ablations were added after the flexible root reached that boundary and are exploratory.

| Model | Status | Test NLL per family |
|---|---|---:|
| asymmetric | scored | 6.86961764761129 |
| baseline | scored | 6.9240962314595595 |
| combined_nb_zero_free | training_failed_numerical_checks |  |
| combined_poisson_zero0 | scored | 6.6735057667478035 |
| hurdle_nb | scored | 6.718433585141261 |
| hurdle_poisson | scored | 6.718433621147111 |
| zero_error0 | scored | 6.936125616474285 |
| zero_error001 | scored | 6.927445503125112 |
| zero_error_free | scored | 6.919038910909843 |
| root01 | scored | 6.718433575838968 |
| root01_asymmetric | scored | 6.647756046868976 |
| root01_zero_error0 | scored | 6.802097808295429 |
| root01_zero_free | scored | 6.605342544347201 |

## Structural and numerical limits

All root-to-tip lengths are 110, up to 8e-9 rounding. A shared homogeneous process and shared observation law therefore imply identical marginal tip distributions. Conditioning on at least one observed tip does not remove this restriction. Real per-species mean counts range from 0.655 to 0.979. No gene identifier is repeated across N11 HOGs. Persistent interspecies discrepancies could reflect lineage-specific evolution, systematic observation differences or ascertainment; these causes are not identified by the count table alone.

See [methods](../../model_improvement_methods.md) for transition equations, root laws, the root-zero identifiability limit, and the distinction between these predictive simulations and a calibrated branch bootstrap. All comparisons use the same family universe; no TE, singleton, large-family or inferred-root filter is applied.


## Interpretation and bootstrap decision

The combined fit substantially improves the bulk occupancy distribution: the
single-species fraction falls from 31.0% under the baseline to 13.9%, compared
with 12.8% observed. Exactly-one-copy profiles increase from 15.1% to 27.4%,
compared with 29.6% observed. Its mean count remains about 0.806 versus 0.754.
The observed fractions lie outside the four descriptive predictive replicates;
these ranges are not calibrated goodness-of-fit p-values. Per-species deviations
remain particularly pronounced. Thus improvement is clear, but adequacy is not
established.

A different starting point with the positive root fixed to exactly one copy
reproduces the combined full-data NLL within 2.5e-8 and closely reproduces its
parameters. This supports a reproducible optimum for that six-category model,
although it does not establish a global optimum or gamma-integration convergence.

The original combined-NB training fit failed convergence after 1,600 iterations;
its root excess mean approached zero. That failure is retained, and its held-out
score is not reported as if it were a converged fit. A simpler exact zero/one-root
training refinement removes the inactive root-tail parameters and is reported
separately as an exploratory candidate. It converged and achieved held-out NLL
per family 6.60534 versus 6.92410 for the baseline (lower is better), a total
held-out log-likelihood improvement of 851.07 across 2,670 test families. This
is predictive improvement, not a calibrated likelihood-ratio test.

The fitted positive-root excess mean approaches zero for both hurdle laws,
including the combined NB model. This is effectively a zero/one-root law; the
NB shape is inactive at that boundary and must not be interpreted. The combined
fit estimates false occurrence at 5.44% and positive-count error at 1.08% in
each direction. These are model-based estimates without independent annotation
validation and may absorb structural misspecification.

Unequal rates help, but the Poisson-root asymmetric fit is start-dependent.
A second start improves NLL from 97290.305 to 96674.225 while root-positive mass,
innovation and error approach zero. Its inclusion probability is about 2.6e-8.
The exact conditional simulator predicts 15.66% single-species and 32.23%
conserved one-copy profiles there. Finite-state cap convergence does not resolve
rare-observation cancellation, boundary identifiability or optimizer multimodality.
This is not evidence that the true innovation rate is zero.

Gamma discretization remains a material numerical approximation. At the optimized
24-category baseline parameters, NLL changes from 97796.913 (24 categories) to
97776.081 (48) and 97767.229 (96), with count-cap checks passing. At the combined
six-category fitted parameters, NLL is 93627.799, 93670.007 and 93685.044 for
6, 12 and 24 categories respectively. These are fixed-parameter integration
checks, not optimized model comparisons. Six-category estimates are exploratory.

The separate innovation/root-zero likelihood slice did not improve toward zero
innovation; it must not be cited as evidence of an empirical flat ridge at the
root01 optimum. The general conditional non-identifiability limit described in
the methods remains a mathematical possibility, illustrated by the separate
asymmetric boundary fit.

**Gate decision:** do not launch a new biological branch bootstrap from these
screening fits. No candidate has established both adequate predictions and stable,
identified parameters. The original biological comparison remains unresolved;
this report makes no new significant expansion/contraction calls. Nominal family
p<0.05 and branch p<0.01/p<0.05 remain the requested reporting thresholds when a
validated branch analysis becomes available. Changing FDR settings cannot repair
these model or numerical issues. The existing bootstrap drivers explicitly reject
the new parameterizations until matching reconstruction/refitting support is
implemented and calibrated.

Next model development should prioritize independent species-level annotation
checks and HOG ascertainment, followed by a small, biologically specified set of
lineage-dependent rates or observation probabilities. A conserved-family component
could address the remaining one-copy peak, but it would not resolve the equal-tip
marginal restriction by itself. These extensions have not been implemented in this
experiment. Avoid fitting a separate unconstrained rate to every branch or choosing
a model because it recovers manuscript HOGs.
