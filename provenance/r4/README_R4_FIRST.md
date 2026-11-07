# ICMPC 2027 R4 — conditional thermal-physics closure and sharp joint-feasibility repair

Status: **working manuscript, not submitted, no materials-process experimental validation.**
Date: 2026-10-30. Original R3 retained separately.

## Validated extensions over R3

1. **Sharp repair:** The exact smallest L1 and L-infinity corrections from the independently reassembled two-node upper/lower interval box corners to any jointly attainable state were derived for all alpha. The minimizer has a closed-form inter-node conductance chosen via a harmonic midpoint. This is a matched-marginal uncertainty-contract result; the general mathematics of interval dependency is not new.
2. **Temperature-dependent physical law:** NIST SRM1450d's certified mean-temperature relation lambda(Tm)=1.10489e-4 Tm is used under a **conditional, unverified local-Fourier constitutive extension** to define the classical Kirchhoff potential (T^2-300^2)/600 K. The original network then solves exactly in transformed potential, while physical absolute temperatures are obtained by square root. Do not call the inherited potential repair distances exact physical-temperature repair distances.
3. **Bounded nominal approximation:** For the hypothetical two-node case, the difference between the transformed potential and true temperature rise is bounded by 50/29403 K (~0.001701 K). This is a mathematical model-form bound, not a measured error, validation, or NIST confidence interval.

## Source legitimacy and limitations

- The only materials-property empirical anchor is NIST SP260-173 (2011), DOI 10.6028/NIST.SP.260-173, thermal conductivity of fibrous-glass board, 280–340 K, 1% expanded uncertainty (k~2). Its effective measured mean-temperature fit is NOT automatically a locally proven Fourier law in the virtual networks.
- Geometry A=0.25 m², L=0.0254 m, path/ambient topology, 20% source bands and even the 1% **fuzzy** coefficient bands are all hypothetical. NIST does not establish fuzzy membership; the use of the same 1% numerical value as a convenient band is not a statistical-to-possibilistic transformation.
- The network's environment is a prescribed 300 K **conducting Dirichlet reservoir**, not a convective Robin boundary; this matters for Kirchhoff linearization.
- No finite-volume or experimental simulation benchmark or universal algorithm-superiority claim is made. Malkawi et al. and other fuzzy papers remain citations for mathematical context; source R11 journal review is independent.
- Kirchhoff transforms are classical: Kim (2001), doi:10.1016/S0735-1933(01)00257-3, and Gama & Pazetto (2021), doi:10.1155/2021/7004581.
- ICMPC 2027 application-type eligibility, detailed 2027 author instructions, paper length and current submission portal require official confirmation. Historic 2026 portals do not apply.

## Reproduce

```bash
python3 r4_exact_sharp_repair.py
python3 reproduce_exact_thermal.py
python3 reproduce_nist_common_k.py
python3 reproduce_matched_marginal_dependence.py
python3 independent_symbolic_proof.py
pdflatex -interaction=nonstopmode -halt-on-error ICMPC2027_R4_Kirchhoff_And_Sharp_Repair_WORKING.tex
pdflatex -interaction=nonstopmode -halt-on-error ICMPC2027_R4_Kirchhoff_And_Sharp_Repair_WORKING.tex
```

Dependencies: Python3 stdlib + sympy, pdflatex; the provided figure PDFs are already vector assets. No web access, field data or restricted files required. Output CSV/JSON files are predetermined checks and the scripts rewrite them with identical results.

The numerical 101-cut grids are **checks**, not a substitute for the analytic all-alpha proofs. The 3,232 R4 flux audits reuse inherited input-corner states; do not add them as distinct experimental cases to R3's 9,696 baseline.
