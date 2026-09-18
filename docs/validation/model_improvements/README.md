# N11 model-improvement experiment: validation and progress

The [plan](../../model_improvement_plan.md) and [methods](../../model_improvement_methods.md)
describe the candidate models and inference limits. All 13,836 N11 families are
retained on the 15-species Vespidae tree. Full-data and training-only comparisons
are currently running; no improved biological significance claim is made yet.

Independent matrix-exponential checks agree with the asymmetric transition kernel
within 7e-15 in tested settings, including extremely small unequal per-copy rates.
Complete small-tree likelihoods agree within 8e-15 for Poisson, hurdle-Poisson and
hurdle-NB roots with separate zero-to-one error. Unbounded simulation frequencies
and supercritical mean checks pass. Root-zero mass cancellation at zero innovation
and zero false occurrence is also verified numerically.

Conditional recovery checks estimate one added parameter at a time while holding
others known. They validate parameter fitting and likelihood replay, not joint
identifiability on the empirical tree. The existing 206 tests / 582 assertions pass.

At identical baseline parameter values, the NLL is 97857.587306 with six gamma
categories, 97879.872895 with twelve and 97852.906639 with twenty-four. This checks
a numerical approximation; each category count needs its own optimized fit before
comparing fitted parameters and predictions.

The predeclared holdout has 11,166 training families and 2,670 test families.
All N11 HOGs belonging to the same original OG remain in one partition; assignments
do not depend on copy counts. Candidate models are defined before held-out scores
are examined. Held-out scores will guide model selection, not provide an unbiased
final estimate of performance after selecting the best candidate on those scores.

The input contains 1,067 single-gene N11 HOGs. Therefore no universal two-gene
minimum is assumed. The exact gene-tree/HOG-construction ascertainment remains
unmodelled. Numerical success alone does not establish empirical model adequacy.
