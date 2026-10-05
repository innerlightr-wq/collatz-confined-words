# Regenerated data

The CSVs for the slope study are not tracked: they total ~25 MB and every one of them is
regenerated deterministically (fixed seeds, exact integer arithmetic) by the scripts in
`shape/code/`. Run them in this order:

```bash
python3 shape/code/stage1_sweep.py        # -> explore_shape.csv
python3 shape/code/stage1b_analysis.py    # -> explore_shape2.csv
python3 shape/code/freeze_predictions.py  # -> ../out/frozen_curve.json
python3 shape/code/stage4_holdout.py      # -> holdout_shape.csv
python3 shape/code/stage8_retest.py       # -> retest_shape.csv
python3 shape/code/stage_b1_explore.py    # -> balance_explore.csv
python3 shape/code/stage_b2_holdout.py    # -> balance_holdout.csv
python3 shape/code/stage_b4_2adic.py      # -> twoadic.csv
```

The derived numbers quoted in the reports are preserved in `shape/out/*.log` and `*.json`, which
*are* tracked, so every claim can be checked without rerunning anything.
