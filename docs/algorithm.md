# Algorithm and scientific scope

A trajectory stores the original sequence and one sequence per substitution. The coding region is an explicit half-open slice `[cds_start, cds_end)`, whose length is a multiple of three. CDS-plus uses one flanking base on each side. The package rejects non-ACGT DNA; data-specific ATG/terminal-stop/annotation filters belong in preprocessing.

## Signature-conditioned events

An SBS96 channel is expanded to its original pyrimidine-centered orientation and reverse complement. A channel weight q contributes b*q to the original orientation and (1-b)*q to its reverse complement. Default b=0.5 is symmetric. For each eligible current sequence position i and alternative nucleotide a, its context supplies weight w(i,a). Event probability is w(i,a) / sum(w). No unchanged-base event or flank target is allowed.

The current sequence is used after every substitution. Hence a mutation can create or destroy a neighboring context, and a site can be hit repeatedly or revert. Profile weights are applied per available opportunity. The resulting 96-channel proportions depend on the sequence and need not equal the profile's original fractions. This is not a model of absolute mutation rate, repair, purifying selection, strand-specific transcription, or mutation timing.

A zero total event weight produces no_available_event; negative/nonfinite weights and malformed profiles are invalid inputs. Packaged signatures are v3.3.1, default GRCh37. Complete canonical SBS96 custom tables can be provided through the library API. The signature loader caches by file path, column, bias, modification time and size; mutation probabilities are rebuilt from the current sequence.

## Neutrality

ProteinIdentity compares original and current translated CDS. The first protein difference terminates the path. For complete paths, K is the number of substitutions before this terminal event, and RUMS is E[K]. K=0 is valid. All prior mutations under this policy are synonymous. Stop-to-stop changes are synonymous under standard-code protein identity; stop-gain/loss changes identity.

BlosumDistance uses global alignment by default. Its distance is original self-alignment score minus original/current alignment score. Distance >= threshold terminates. The baseline is always the original sequence; distances are not summed per mutation. Intermediate protein changes can be neutral under this criterion, so the result is criterion-neutral stopping time. A positive BLOSUM threshold is a modeling parameter, not a validated experimental functional cutoff. The explicit stop-change policy can force termination regardless of distance.

A custom callable sees original CDS, current CDS and the mutation event. It returns bool or EvaluationResult. Incomplete/fixed-step terminal reason names are reserved. No-event or max_steps paths are not complete stopping-time samples. The scalar stopping-time aggregate refuses them rather than censoring or discarding silently; explicit survival analysis is a separate downstream task.

Fixed-n endpoint neutrality is a separate metric: phenotype equality after exactly n steps, including possible restoration after intermediate non-neutral states. It is not the survival probability of never leaving neutrality. `by_score=True` retains the older reciprocal mean-distance exploratory measure and returns infinity at mean distance zero; it is not a neutral probability.

## Reproducibility and limits

Generator state is local to each trajectory; SeedSequence creates independent indexed streams for batched trajectories. CLI streams include stable sequence ID and signature ID in their seed derivation. Inputs, configuration, package and dependency versions should be archived with each scientific run. No model results establish genome adaptation, longevity causation or tumor fitness without suitable biological/null-model analyses.
