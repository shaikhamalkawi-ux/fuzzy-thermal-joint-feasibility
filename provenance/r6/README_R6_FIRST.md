# ICMPC 2027 R6 — Exact Three-Node All-α Joint Repair

**31 October 2026 — 9-page research working manuscript, not submitted.**

## Source lock

The 8-page ICMPC R5 NIST15 package is preserved unmodified as `UPSTREAM_R5_IMMUTABLE.zip`; it remains the scientific baseline for the two-node, matched-conductance, NIST certified-fit/historical-data, and conditional Kirchhoff results. No additional materials measurements, geometry calibration, process data, fuzzy-membership calibration, or ICMPC venue permission have been obtained.

## New scientific contribution

For the original synthetic three-node chain with *fixed* grounded ambient conductances `g_i=1 W/K`, heat source cuts `(5±q)*(2,1,2/5) W`, and independent inter-node conductances `a,b∈[1−q/5,1+q/5]`, exact nodal-balance and whole-domain polynomial inequalities establish the **attained** minimum corrections to upper/lower marginal-envelope corners under L1 and L∞ at **every** α. The previous R5 supplied only three-node envelope gaps, not their nearest jointly valid states.

Let G_+(q)=sum_i U_i−(17/5)(5+q), G_-(q)=(17/5)(5−q)−sum_i L_i, both in watts. Then d1_±(q)=G_±/(1 W/K) and d∞_±(q)=G_±/(3 W/K), in K. The equality in L∞ is **not automatic** for general thermal graphs; it is established by an explicit 3-node admissible equal-repair construction, with 8 exact coefficient-domain inequalities and nonzero positive integer lower bounds on all normalized polynomials. All 202 exact Fraction cases are independent consistency checks, not the continuum proof.

At q=1:

| Extremum | d1 (K) | d∞ (K) |
|---|---|---|
| Upper | 1855736/2541891 ≈0.730061 | 1855736/7625673 ≈0.243354 |
| Lower | 3695632/7625673 ≈0.484630 | 3695632/22877019 ≈0.161543 |

These are mathematical certificates for a declared illustrative network; no measured material heat balance is claimed. The known joint-versus-marginal dependence problem and elementary conservation/norm inequalities are not claimed as novel universal theorems.

## Reproduce

```
python3 R6_3node_sharp_repair_certificate.py
python3 R6_independent_fraction_replay.py
python3 build_R6_paper.py
pdflatex -interaction=nonstopmode -halt-on-error ICMPC2027_R6_ThreeNode_SharpRepair_SourceAudited_WORKING.tex
pdflatex -interaction=nonstopmode -halt-on-error ICMPC2027_R6_ThreeNode_SharpRepair_SourceAudited_WORKING.tex
```

Dependency: Python `sympy` and `matplotlib` for the first verifier and plot; the second verifier uses only Python standard library. The `.tex` reuses PDF figures included in this package.

## Submission gate

Do not submit until the 2027 ICMPC paper type, exact submission route and final template are verified; and materials-process physical scenario/geometry/load relevance and direct comparable prior art are established. Preserve the separate 99fuzzy R12 mathematical review, and avoid dual-publication of overlapping examples.
