#!/usr/bin/env python3
from pathlib import Path
import csv, json, shutil, re, textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
R4=ROOT.parent/'r4'
# Carry original proven R4 source and vector figures unchanged.
for name in ('Figure_2node_Energy_Balance_Gap.pdf','Figure_NIST_MatchedMarginal_Dependence.pdf'):
 shutil.copyfile(R4/name,ROOT/name)
source=(R4/'ICMPC2027_R4_Kirchhoff_And_Sharp_Repair_WORKING.tex').read_text(encoding='utf8')
cert=json.loads((ROOT/'NIST15_source_audit_certificate.json').read_text())
with (ROOT/'NIST_1450d_TableI1_residuals_recomputed.csv').open(newline='') as f:
 rows=list(csv.DictReader(f))
# Manuscript-quality vector provenance plot -- visualization of NIST source residuals, NOT a fitted model.
fig,ax=plt.subplots(figsize=(6.3,2.3))
xs=[float(r['Tmean_K']) for r in rows]
ys=[float(r['relative_residual_pct']) for r in rows]
ax.scatter(xs,ys, s=38, marker='o')
ax.axhline(1.0,linestyle='--',linewidth=0.9)
ax.axhline(-1.0,linestyle='--',linewidth=0.9)
ax.axhline(0,linewidth=0.65,alpha=0.55)
ax.set(xlabel='Published specimen mean temperature (K)', ylabel='Relative source residual (%)',
       ylim=(-1.18,1.18), xlim=(272,348))
ax.grid(axis='y',alpha=0.15)
fig.tight_layout()
fig.savefig(ROOT/'Figure_NIST_2014_TableI1_SourceResiduals.pdf')
plt.close(fig)
source=source.replace('When Marginal Thermal Bounds Lose Joint Feasibility: Exact Repair and Material-Property Dependence Certificates',
    'When Marginal Thermal Bounds Lose Joint Feasibility: Exact Repairs and a Source-Audited Material Sensitivity')
source=source.replace('30 October 2026 -- ICMPC 2027 R4 scientific working draft; not submitted',
    '31 October 2026 -- ICMPC 2027 R5 scientific working draft; not submitted')
source=source.replace('The exact certificates quantify representational dependence, not empirical processing accuracy, solver superiority, or priority over established interval and Kirchhoff methods.',
    'A separate audit of 15 published SRM 1450d measurements shows finite-sample relative deviations from the certified mean-temperature fit below 0.668\\%, and compares a descriptive observed-residual sensitivity scenario against the original illustrative 1\\% band. These data do not validate the network or specify fuzzy membership. The certificates quantify representational dependence rather than empirical process accuracy or solver superiority.')
insert=r'''
\subsection{Source-level check using the NIST 1450d measurement table}
\label{sec:sourceaudit}
The certified relation in Eq.~\eqref{eq:nist} anchors the conductivity scale, but its empirical source contains more information than the single coefficient. Zarr, Heckert and Leigh's later retrospective paper~\cite{nist2014data} republishes the 15 source measurements for SRM 1450d Lot 2009 (Table I1), with individual conductivity, density, mean-temperature, plate temperature difference and thickness. We independently recomputed the residuals
\begin{equation}\label{eq:source_resid}
 e_j=\lambda_{\mathrm{exp},j}/(1.10489\times10^{-4} T_{m,j})-1,\qquad j=1,\ldots,15.
\end{equation}
There are three specimens at each nominal mean-temperature group 280, 295, 310, 325 and 340~K. The published 15 decimal measurements give mean absolute $|e_j|$ of 0.27063\%, relative RMSE $\sqrt{15^{-1}\sum_j e_j^2}$ of 0.34034\%, and $\max_j|e_j|=0.66776\%$ (row 15 at printed $T_m=339.999$~K). All 15 absolute deviations are below 1\% numerically. Two printed mean temperatures are 279.998~K, 0.002~K below the nominal 280~K lower certification boundary; they are retained and flagged rather than rounded without disclosure. These historical source observations \emph{are not a held-out validation sample}, and the certified 1\% expanded uncertainty ($k\approx2$) is neither a deterministic guarantee on every measurement nor a fuzzy membership function.
\begin{table}[t]
\caption{Published NIST SRM 1450d Lot-2009 measurements, recomputed against its certified mean-temperature relation. Three measurements per nominal temperature group. Groupwise values are descriptive, not model validation or uncertainty quantiles.}\label{tab:nist_rows}
\centering\small
\begin{tabular}{@{}rccc@{}}\toprule
Nominal $T_m$ (K) & $n$ & Measured $\lambda$ range (W/m K) & Max. $|e_j|$ (\%) \\\midrule
280 & 3 & 0.030898--0.031123 & 0.6022\\
295 & 3 & 0.032450--0.032614 & 0.4429\\
310 & 3 & 0.034089--0.034283 & 0.4750\\
325 & 3 & 0.035723--0.035988 & 0.5175\\
340 & 3 & 0.037524--0.037817 & 0.6678\\\bottomrule
\end{tabular}
\end{table}
\begin{figure}[t]
\centering\includegraphics[width=.74\linewidth]{Figure_NIST_2014_TableI1_SourceResiduals.pdf}
\caption{Residuals of 15 historic NIST thermal-conductivity measurements relative to the original 1450d certified mean-temperature fit. Dashed $\pm1\%$ lines denote the published expanded-uncertainty percentage for context only; they are not a fuzzy cut or a guaranteed per-specimen error bound.}\label{fig:nist_rows}
\end{figure}

For a \emph{descriptive sensitivity}, set $\eta_{\mathrm{obs}}=\max_j|e_j|=0.0066775672318$ in the already derived matched-marginal model of Eq.~\eqref{eq:matched}, keeping the source interval half-width $\sigma=0.2$ purely illustrative. The observed-residual percentage is \emph{not} fitted as a probability/possibility distribution, is not guaranteed to enclose unseen samples, and does not imply independently uncertain paths in the tested NIST boards. The calculation only probes sensitivity of the hypothetical network to a smaller relative width:
\begin{center}\small
\begin{tabular}{@{}lrr@{}}\toprule
Assumed $\eta$ at $\alpha=0$ & Apparent upper excess (W) & Apparent lower deficit (W)\\\midrule
$\eta=0$ (crisp control) & 0 & 0\\
$\eta_{\mathrm{obs}}=0.00667757$ (descriptive) & 0.000289826 & 0.000194080\\
$\eta=0.01000000$ (R4 illustration) & 0.000433551 & 0.000290967\\\bottomrule
\end{tabular}
\end{center}
The corresponding attained sharp $L_1$ correction distances for the descriptive scenario are 0.00089434~K (upper) and 0.00059094~K (lower), within the \emph{potential-temperature} scale of the declared model. The positivity of the artificial-box gaps under independent paths and their exact disappearance under a shared path factor are unchanged by the choice of $0<\eta\le0.01$: this is an analytic dependence result, not an inference from the 15 observations. The complete 15-row transcription, published source locators and exact-fraction calculation are supplied separately to make this limited evidence step auditable. The original NIST raw-data ZIP is listed by NIST but is not redistributed or claimed acquired here.
'''
needle=r'\section{Sharp repair of the independently recombined extreme state}'
assert needle in source
source=source.replace(needle,insert+'\n'+needle)
# Key limitation paragraph updated with added within-source evidence
source=source.replace('This is \\textbf{not} a measurement-based materials-processing case',
     'This is \\textbf{not} a measurement-validated materials-processing case')
source=source.replace('In the control, only nominal conductivity is anchored to NIST; the geometry, load boxes, and fuzzy-support mapping are synthetic and have no uncertainty calibration.',
    'In the control, NIST provides a certified conductivity formula and 15 historical measurements for a separate finite-sample consistency audit; the network geometry, source boxes, and fuzzy-support mapping are still synthetic and have no uncertainty calibration.')
source=source.replace('A conditional NIST-linked Kirchhoff potential preserves the distinction under a temperature-dependent Fourier law and quantifies the fixed-conductivity approximation.',
    'A conditional NIST-linked Kirchhoff potential preserves the distinction under a temperature-dependent Fourier law. A 15-measurement source audit characterizes the published conductivity-fit residuals and tests a narrower descriptive sensitivity without asserting network validation.')
source=source.replace('\\begin{thebibliography}{9}', '\\begin{thebibliography}{12}')
needle_bib='\\bibitem{kim2001}'
new_bib=r'''\bibitem{nist2014data} R.~R. Zarr, N.~A. Heckert, and S.~D. Leigh, ``Retrospective Analysis of NIST Standard Reference Material 1450, Fibrous Glass Board, for Thermal Insulation Measurements,'' \emph{J. Res. Natl. Inst. Stand. Technol.}, vol.~119, pp.~296--370, 2014, Appendix I, Table I1. \href{https://doi.org/10.6028/jres.119.012}{doi:10.6028/jres.119.012}.
'''
source=source.replace(needle_bib,new_bib+needle_bib)
assert source.count('\\bibitem{nist2014data}')==1
out=ROOT/'ICMPC2027_R5_NIST15_SourceAudited_WORKING.tex'
source=source.replace('The grid supports the mathematics but does not establish its continuum validity, which follows from the displayed monotonicity/derivative argument.','The grid supports the mathematics but does not establish its continuum validity, which follows from the displayed monotonicity/derivative argument. The additional R5 NIST-Table-I1 audit script recomputes all 15 conductivity-fit residuals using exact decimal fractions; a separate exact checker verifies 3,232 mathematical source/conductance corner states over 101 alpha levels for the descriptive observed-residual scenario. These are deliberately not treated as new physical observations or as an independent experimental test.')
out.write_text(source,encoding='utf8')
print('wrote',out, 'chars',len(source))
