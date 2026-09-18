# N11 model-improvement experiment

Status: screening and numerical checks completed; final comparison records the
held-out and second-start results and the bootstrap gate decision. Primary data: all 13,836 observed N11 families on the matching
15-species Vespidae dated tree. N0/17-species results remain an archived benchmark.
No count-size, singleton or inferred-root filter will be applied.

## Ordered experiments

1. Check the six-category gamma approximation against 12 and 24 categories at
   identical parameters, then re-fit if the likelihood changes materially.
   Keep count-cap checks separate from gamma integration checks.
2. Implement and independently validate separate duplication and loss rates,
   retaining a single global innovation rate. Gamma scales both per-copy rates.
   Validate transition probabilities against a large-state matrix exponential,
   critical/pure-death/pure-immigration limits and simulation moments.
3. Separate zero-to-one annotation error from positive-count error. Re-fit with
   zero error fixed to zero, 0.001 and separately estimated; these are sensitivity
   assumptions, not claims that an externally measured error rate is available.
4. Compare the Poisson root with hurdle shifted-Poisson and hurdle shifted-NB
   roots. Estimate root-zero mass independently of positive ancestral counts.
   Include the exact zero/one-root boundary (zero excess mean) so an optimizer
   drifting toward that limit is not mistaken for a well-identified positive tail.
5. Compare changes individually and jointly. Screen with multiple starts and
   whole-dataset predictive replicates. Require convergence, enlarged-cap checks
   and gamma-category stability. Use a reproducible family holdout for finalist
   predictive likelihood comparison, grouping N11 HOGs by their original OG identity so related subfamilies remain
   in one partition. Identical count profiles alone do not imply dependence.
6. Audit the supplied HOG catalogue and ascertainment. Document known inclusion
   rules and unresolved mechanisms. Do not invent a correction from observed
   occupancy or remove families to manufacture agreement. Any implemented
   selection rule must be identical in likelihood and simulation.
7. For an adequate, identifiable finalist, run simulation recovery and refitted
   bootstrap calibration before biological significance. Retain nominal family
   p<0.05 and branch p<0.01/0.05 reporting. Do not use marginal-MAP transition
   surprise alone as a directional test. If no candidate is adequate, report
   that outcome and its remaining mismatch rather than publish new branch calls.

## Model comparison targets

Assess species occupancy, exactly-one-copy-in-every-species fraction, singleton
and single-gene fractions, mean counts, count-difference tails, per-species means
and zero fractions. Prefer held-out predictive performance and model adequacy;
manuscript HOG recovery is not a model-selection objective. AIC comparisons, if
reported, must use identical observations and conditioning. Boundaries, weakly
identified parameters and optimizer failures must be explicit.

Screening uses small numbers of full-dataset simulations; these are predictive
checks, not calibrated branch p-values. Longer bootstrap work is conditional on
passing the adequacy gate. All scripts, seeds, fit logs and comparison tables
will be retained. The current working branch is innovation-bdi.
