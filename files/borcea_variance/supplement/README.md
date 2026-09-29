# Exact finite verification

The accompanying `../borcea_variance.pdf`, especially Appendices A and B, states the mathematical inequalities checked here.

```bash
python3 self_test.py
python3 verify_all.py --bits 100
python3 verify_all.py --bits 160
```

These commands require Python 3 and its standard library only. Use `python` on systems where that names Python 3. The wrapper works from any current directory because it resolves its input files relative to its own location. Reports are written into `rerun_logs/`.

| Certificate | Degrees | Blocks | Leaves |
|---|---|---:|---:|
| `certificate_n4_to_n13.json` | 4–13 | 10 | 12,821 |
| `certificate_n14_to_n100000.json` | 14–100000 | 62 | 20,850 |

`verify_low.py` uses exact arithmetic utilities from `verify_mid.py`; these are not fully independent checker implementations. `verify_all.py` runs both and checks their degree endpoints. The two complete runs at 100 and 160 bits are included as JSON logs. Running at both precisions is a useful cross-check, not a formalization of the mathematical lemmas.

Optional generators are `generate_low.py`, `propose_mid.py`, `regenerate_low.py`, and `regenerate_mid.py`; their command-line help explains usage. They require NumPy and can use floating-point proposals. Their decisions are not trusted for acceptance: the exact checkers read only the resulting partition and rule identifiers and recompute each accepted inequality. The present revision reran acceptance checks, not the entire proposal generation process.

The source root's `audit/independent_checks.py` separately checks all finite covers by sorted binary addresses and rational volume sums. It does not import these modules and does not independently reimplement every numerical inequality.

The mathematical proof for degrees at least 100001 is Proposition 4.1 of the paper; no finite run is extrapolated to that range. The ten small-degree leaves that use the analytic near-equality argument rely on Lemma 6.6, proved in the paper.

For reproducibility, Python 3.10 or later is recommended. Run the diagnostic commands without `-O`, so their assertion checks remain enabled.
