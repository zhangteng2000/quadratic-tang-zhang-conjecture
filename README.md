# Quadratic Tang--Zhang Conjecture

This repository contains the manuscript and the complete computational
materials accompanying a computer-assisted proof of the quadratic
Tang--Zhang conjecture.

## Main result

Let \(p\) be a complex polynomial of degree \(n\ge 2\), all of whose zeros
lie in the closed unit disk. Let \(\zeta\) be a zero of \(p\), and let
\(c_1,\ldots,c_{n-1}\) be the critical points of \(p\), counted with
multiplicity. We prove that

\[
\sum_{j=1}^{n-1}\frac{1}{|\zeta-c_j|^2}\ge n-1.
\]

The equality case is also determined.

As a consequence of the power-mean inequality, the corresponding estimate
holds for every exponent \(\lambda\ge 2\).

## Manuscript

- `quadratic_tang_zhang_revised.pdf` — compiled manuscript.
- `quadratic_tang_zhang_revised.tex` — complete LaTeX source.

## Computational supplement

The directory `supplement/` contains the computer-assisted part of the proof.

In particular:

- `certificate.json` — rigorous interval certificate for the finite range.
- `certify.py` — certificate generator.
- `verify_certificate.py` — certificate verifier.
- `verify_independent.py` — separately implemented certificate verifier.
- `analytic_audit.py` — exact audit of the numerical inequalities used in
  the analytic part of the proof.
- `verification_original.json` — recorded output of the first verifier.
- `verification_independent.json` — recorded output of the second verifier.
- `verification_analytic.json` — recorded output of the analytic audit.

The file `SHA256SUMS.txt` records SHA-256 checksums for the distributed files.

## Reproducing the verification

Only Python 3 and the Python standard library are required.

From the repository root, run

```bash
cd supplement
python verify_certificate.py certificate.json
python verify_independent.py certificate.json
python analytic_audit.py
```

The supplied certificate covers the finite parameter range used in the paper.
The large-degree range is treated analytically in the manuscript.

## Formalization status

A Lean 4 formalization of the proof is in preparation.

The Lean formalization is not yet part of the present release. It will be
added to this repository when the corresponding source files are ready for
public verification.

## Verification status

The supplied computational certificates have been replayed successfully by
the verification programs included in this repository.

The repository is made public to facilitate independent checking of both the
mathematical argument and the computer-assisted component.

## Additional files

- `revision_report_zh.md` — detailed record of the mathematical and editorial
  revisions.
- `reference_audit_zh.md` — bibliographic, DOI, and source-location audit.
- `editorial_audit.json` — automated checks of cross-references and
  typesetting consistency.

These files document the preparation and verification of the manuscript and
are not required for running the proof certificate.

## Citation

When citing this work, please cite the manuscript and, when appropriate, the
specific GitHub release used for verification.

A permanent archival DOI will be added once a fixed release has been deposited
in an archival repository.

## Author

Teng Zhang  
School of Mathematics and Statistics  
Xi'an Jiaotong University  
Xi'an 710049, P. R. China

## Lean 4 formalization

A complete Lean 4 formalization of the main theorem and equality
classification is available at:

https://github.com/zhangteng2000/quadratic-tang-zhang-lean4
