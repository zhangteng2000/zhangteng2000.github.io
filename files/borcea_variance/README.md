# Borcea's 2-variance conjecture

Teng Zhang. Manuscript and supplementary materials, 29 September 2026.

- [Paper (PDF)](borcea_variance.pdf)
- [LaTeX source](borcea_variance.tex)
- [Code, certificates and verification reports](supplement/)
- [Complete download](../borcea_variance_2026-09-29.zip)
- [Mathematical and computational audit](AUDIT_REPORT.md)
- [Bibliographical audit](BIBLIOGRAPHY_AUDIT.md)
- [100-bit run log](audit/delivered_100bit.txt) and [160-bit run log](audit/delivered_160bit.txt)

## Reproduce the finite verification

Use Python 3.10 or later. From this directory:

```bash
cd supplement
python3 self_test.py
python3 verify_all.py --bits 100
python3 verify_all.py --bits 160
```

On Windows, use `python` if that is the name of your Python 3 command. The acceptance checks use only the Python standard library. Optional certificate generators require NumPy. Reports are written to `supplement/rerun_logs/`. Run without `-O` so that diagnostic assertions remain enabled.

| Certificate | Degrees | Degree blocks | Leaf boxes |
| --- | --- | ---: | ---: |
| `certificate_n4_to_n13.json` | 4-13 | 10 | 12,821 |
| `certificate_n14_to_n100000.json` | 14-100000 | 62 | 20,850 |
| Total | 4-100000 | 72 | 33,671 |

The paper supplies the mathematical reductions and the analytic argument for degrees at least 100001. The two precision runs use the same checking algorithms; they do not constitute two independent proof implementations or a proof-assistant formalization. The separate coverage diagnostic is `audit/independent_checks.py`.

## Citation

T. Zhang, *Supplementary code and certificates for Borcea's 2-variance conjecture*, version 2026-09-29, 2026.

See [CITATION.bib](CITATION.bib) for BibTeX. The manuscript cites this archive as [Zha26S]. `SHA256SUMS.txt` records the exact file hashes for this version.

## Compile the manuscript

```bash
pdflatex -interaction=nonstopmode -halt-on-error borcea_variance.tex
pdflatex -interaction=nonstopmode -halt-on-error borcea_variance.tex
```

The manuscript has an embedded bibliography; BibTeX is not required. The public version adds the code citation and the author-supplied acknowledgments and AI tools disclosure at the end of the Introduction. The verification source files and certificates are unchanged.

The audit reports and delivered run logs document the supplied manuscript before the public code citation and the author-supplied disclosure were added. The upload preparation also reran the complete 100-bit finite verification successfully; those reports are in `audit/preupload_100bit/`.
