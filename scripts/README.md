# Scripts

Gene-symbol drug screen. No genome.

```bash
python scripts/screen.py --out scripts/out
```

Asks Open Targets and ChEMBL what approved or late-phase drugs are
known for *BUB1B* and neighbouring spindle-assembly-checkpoint genes,
then pulls openFDA label flags for unique approved names.

`hits.tsv` is the working table. `class` is a first-pass harm filter.
Ranking and paediatric judgement stay in the Track 2 report.
