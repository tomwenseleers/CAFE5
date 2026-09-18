# Illustrative count-table benchmark

`counts.tsv` and `tree.nwk` are a fixed example from public Vespidae HOG data.
See `provenance.json` for sources, selection and hashes. These files demonstrate
software behavior and likelihood comparison; they are not a curated biological
reference or a claim that every retained family lacks TE association.

From the repository root:

```bash
python3 scripts/innovation/analyze.py \
  examples/innovation_benchmark/counts.tsv \
  examples/innovation_benchmark/tree.nwk --fit-only -o example_fit
```

The full example is computationally substantial. The small generated example in
`scripts/innovation/validate_workflow.py` is suitable for a quick installation test.
