# Quadratic Tang–Zhang manuscript and reproducibility package

This package contains the revised English research manuscript, its compiled PDF,
three-round internal review notes, a reference-verification ledger, and the
computational supplement used in the proof.

## Main files

- `quadratic_tang_zhang_revised.tex`: complete editable English LaTeX source.
- `quadratic_tang_zhang_revised.pdf`: compiled 17-page manuscript.
- `revision_report_zh.md`: mathematical and editorial review, in Chinese.
- `reference_audit_zh.md`: source locations, DOIs, and access limitations.
- `editorial_audit.json`: automatic cross-reference and typesetting checks.
- `supplement/`: certificate, generator, two checkers, and recorded outputs.

## Compile

A TeX installation containing the packages in the supplied amsart preamble is
required. In particular, bbm, TikZ/PGFPlots, doi, bookmark, and microtype must be
installed. No external images or bibliography database are required.

```sh
pdflatex -interaction=nonstopmode -halt-on-error quadratic_tang_zhang_revised.tex
pdflatex -interaction=nonstopmode -halt-on-error quadratic_tang_zhang_revised.tex
```

## Reproduce the finite verification

Use Python 3 with its standard library. Run from the supplement directory:

```sh
cd supplement
python verify_certificate.py certificate.json
python verify_independent.py certificate.json
python analytic_audit.py
```

The certificate covers every integer `5 <= m <= 999999` and the continuous
parameter interval `0 <= a <= 1`. The actual endpoints `a=0,1` are treated
analytically in the manuscript. The theorem uses `n=m+1`.

The original checker imports its arithmetic implementation from `certify.py`.
The separate checker imports neither original program and uses exact rational
coefficients and 160-bit directed power enclosures. Both verify the same
mathematical reduction. Their agreement is not a Lean formalization or independent
human peer review.

The original checker and generator retain their supplied implementation; no
mathematical acceptance tests were loosened to obtain the equality statement.
The extension from `s<1` to `s<=1` is justified in the manuscript, not by changing
certificate thresholds.

## Regenerate, rather than only replay

```sh
python certify.py --stop 1000000 --out regenerated_certificate.json
python verify_independent.py regenerated_certificate.json
```

Regeneration need not produce a byte-identical search history. An accepted
certificate must pass all inequality and complete-cover tests. Ordinary floating
point in the generator is used only for proposals or reporting, never for the
acceptance of a rectangle.

## Scope and publication status

This is an AI-assisted research manuscript. It is not a published theorem,
external referee report, or proof-assistant-certified result. Its claims are the
quadratic inequality, equality classification, and the power-mean consequences
for exponents at least two. It does not claim the exponent-one conjecture.

Before submission/publication, the author should arrange independent review of
the analytic reduction and a permanent public archive for the essential
computational supplement. No repository DOI has been invented for this package.
