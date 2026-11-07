# ICMPC 2027 R5 — Source-audited material sensitivity (31 October 2026)

**Latest manuscript:** `ICMPC2027_R5_NIST15_SourceAudited_WORKING.pdf` (8-page working draft; not submitted).
**Previous baseline:** R4 ZIP is embedded byte-for-byte as `UPSTREAM_R4_IMMUTABLE.zip`. No R4 mathematical conclusions are altered.

## New empirical material-property evidence (not model validation)

- NIST JRES 119 (2014), DOI `10.6028/jres.119.012`, printed page 360, Appendix I, Table I1: 15 measurements from SRM 1450d Lot 2009.
- NIST SP260-173 (2011), DOI `10.6028/NIST.SP.260-173`: lambda_ref(Tmean)=0.000110489*Tmean W/(m K), certified for 280–340 K mean temperature, expanded measurement uncertainty 1% (k~2).
- `NIST_1450d_TableI1_source_rows.csv`: 15 manually transcribed published table rows; the official downloadable dataset ZIP was not successfully obtained and **is not claimed archived**.
- `NIST_1450d_TableI1_residuals_recomputed.csv`: all 15 source relative residuals and exact rational representation.
- `NIST15_source_audit_certificate.json`: descriptive statistics and scope guardrails. Max observed abs residual 0.6677567%, relative RMSE 0.3403380%, all 15 deviations below 1%.
- `NIST_sensitivity_results.csv`: sensitivity of the unchanged matched-marginal mathematical model to eta=0, max observed descriptive residual (0.006677567), original illustrative eta=0.01.
- `independent_3232_cut_checker.py`: 3,232 exact-Fraction hypothetical thermal states at 101 alpha levels for descriptive-eta scenario. All actual states conserve energy; gap applies only to unattainable marginal recombinations.

## Reproduce

```
python3 reproduce_nist15_source_sensitivity.py
python3 independent_3232_cut_checker.py
pdflatex -interaction=nonstopmode -halt-on-error ICMPC2027_R5_NIST15_SourceAudited_WORKING.tex
pdflatex -interaction=nonstopmode -halt-on-error ICMPC2027_R5_NIST15_SourceAudited_WORKING.tex
```

R4 full upstream methods/scripts are embedded. This file is source data and a computed audit, **not** an independent materials-processing experiment or a fuzzy uncertainty calibration. NIST's 1% expanded uncertainty is not a hard finite-sample or interval guarantee. The entire network geometry, heat inputs, fuzzy membership widths and local-Fourier closure are conditional and synthetic.

**Submission HOLD**: materials-process case/geometry/load evidence, source-matched thermal-UQ comparator, event official 2027 submission format and article-type acceptance, and author approvals.
