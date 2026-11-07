# Fuzzy thermal networks: exact joint feasibility and attained repairs

**Research-status notice:** Source-defined, synthetic two-/three-node steady-state conduction examples, not a field-calibrated materials-processing predictor, a solver performance comparison, or an industrial safety/control certificate. Scientific model baseline **R6 (2026-10-31)**. A separate editorial manuscript candidate **R6P12 (2026-11-08)** remains unpublished/unsubmitted.

## Mathematical question

Individual exact temperature extrema can be attained with different admissible heat-input and conductance configurations, but their simultaneous recombination need not belong to the jointly attainable state set. The existing project derives precise conservation gaps and sharp attained minimum L1 and L-infinity repairs for stated networks, including an all-alpha three-node link-admissibility certificate. The general joint-set-versus-box distinction, interval hulls, monotonicity approach and Kirchhoff transform are established prior art, not claimed as new universal methods.

## Quick reproduction

```bash
python -m pip install -r requirements.txt
python run_verify.py
```

`run_verify.py` uses the original byte-identical **released code files** in `scripts/`, runs them in clean temporary directories with source-defined inputs in `data/`, and compares newly generated CSV/JSON results to the original R4/R5/R6/G2 outputs. It does not require network access, proprietary data, manual coders, or a physical experiment. Python 3 and `sympy==1.14.0` are required. For optional LaTeX PDF compilation: `pdflatex` and the figure assets in `manuscript/` are needed. The repository's reproducibility test contains 19 expected output-parity checks: R4=9, R5=5, R6=3, G2=2. The continuum certificate depends on exact symbolic identities and positivity tests, not on the grid checks alone.

## Repository layout

- `scripts/r4`, `scripts/r5`, `scripts/r6`, `scripts/g2`: frozen Python reproduction and independent-certificate scripts, readable separately.
- `data/r4`, `data/r5`, `data/r6`, `data/g2`: numerical tables, source-limited transcriptions, exact outputs and rational-certificate artifacts.
- `figures/`: four source-defined author-generated plot assets (not external publisher PDFs). The editable manuscript/working-paper PDFs are held separately for author and venue policy approval, and are not in the default public repository staging. **No accepted conference article is claimed.**
- `provenance/ORIGINAL_FILE_SHA256.json`: file-by-file source lineage, original archive attribution and hashes.
- `run_verify.py`: clean replay from published source files. A successful run checks byte equality of 19 outputs; it does not establish empirical validation.

## Data provenance and restrictions

The `NIST_1450d_TableI1_source_rows.csv` is a transcription of **15 historical published NIST measurements** from Zarr, Heckert and Leigh (2014), Appendix I, Table I1 (https://doi.org/10.6028/jres.119.012). NIST SP 260-173 (https://doi.org/10.6028/NIST.SP.260-173) supplies the conductivity formula. The values are used only as source checks and descriptive sensitivity—not as a calibrated fuzzy membership distribution or independent network validation. The original NIST PDFs, publisher full texts and paywalled research papers are **not redistributed** here.

Every source-dependent claim has a model/data boundary. The project uses no private industrial temperatures, no physical field geometries, and no human coding or participant data.

## Citation, license, and manuscript status

This is a **public inspection/reproducibility repository**. The public source repository is https://github.com/shaikhamalkawi-ux/fuzzy-thermal-joint-feasibility. Use a specific Git commit SHA for reproducibility. No software DOI is claimed. A `CITATION.cff` will be added **after actual software authorship is confirmed**; no invented authors or DOI have been assigned. No open-source license has been granted yet: public visibility does not itself grant permission to reuse or redistribute the material. The rights status and conference self-archiving requirements must be checked before publishing any manuscript or its editable TeX. Working manuscript files are intentionally excluded from the public-code staging bundle. See `RIGHTS_AND_RELEASE_STATUS.md` and `CITATION_GUIDANCE.md`.

Contact and institutional information are intentionally not invented. No submission, publisher endorsement, or proceedings indexing is asserted.
