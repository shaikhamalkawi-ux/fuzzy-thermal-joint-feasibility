#!/usr/bin/env python3
"""Editorial build of ICMPC R6 from unchanged R5 LaTeX and new exact R6 proof."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42

root=Path(__file__).resolve().parent
old=(root/'R5_upstream/ICMPC2027_R5_NIST15_SourceAudited_WORKING.tex').read_text(encoding='utf-8')
assert old.count('\\section{Conditional temperature-dependent conductivity closure}')==1

# Figure uses exactly reconstructed upper and lower L1, Linf distances.
with (root/'R6_3node_101cut_sharp_repair.csv').open(encoding='utf-8') as f:
 rows=list(csv.DictReader(f))

def value(x):
 p,n=x.split('/');return int(p)/int(n)
up=sorted((r for r in rows if r['corner']=='upper'),key=lambda r:float(value(r['alpha'])))
lo=sorted((r for r in rows if r['corner']=='lower'),key=lambda r:float(value(r['alpha'])))
a=[value(r['alpha']) for r in up]
fig,ax=plt.subplots(figsize=(6.9,3.6))
ax.plot(a,[value(r['d1_K']) for r in up],label=r'Upper $L_1$ distance',linewidth=1.9)
ax.plot(a,[value(r['d1_K']) for r in lo],label=r'Lower $L_1$ distance',linewidth=1.9)
ax.plot(a,[value(r['dinf_K']) for r in up],linestyle='--',label=r'Upper $L_\infty$ distance',linewidth=1.7)
ax.plot(a,[value(r['dinf_K']) for r in lo],linestyle='--',label=r'Lower $L_\infty$ distance',linewidth=1.7)
ax.set(xlabel=r'Membership level $\alpha$',ylabel='Minimum jointly feasible correction (K)',xlim=(0,1),ylim=(0,None))
ax.grid(alpha=.2)
ax.legend(frameon=False,fontsize=8,ncol=2)
fig.tight_layout()
fig.savefig(root/'Figure_R6_ThreeNode_SharpRepair.pdf',bbox_inches='tight')
fig.savefig(root/'Figure_R6_ThreeNode_SharpRepair.png',dpi=165,bbox_inches='tight')
plt.close(fig)

oldtitle='When Marginal Thermal Bounds Lose Joint Feasibility: Exact Repairs and a Source-Audited Material Sensitivity'
newtitle='When Marginal Thermal Bounds Lose Joint Feasibility: Exact Two- and Three-Node Repair Certificates and Material Sensitivity'
assert oldtitle in old
old=old.replace(oldtitle,newtitle,1)
old=old.replace('ICMPC 2027 R5 scientific working draft','ICMPC 2027 R6 scientific working draft',1)
old=old.replace('For the independent-path case, we derive \\emph{sharp} all-level minimum $L_1$ and $L_\\infty$ repairs to both unattainable extreme vectors and construct the conductance that attains them.',
'For the matched-marginal independent-path case, we derive \\emph{sharp} all-level minimum $L_1$ and $L_\\infty$ repairs to both unattainable two-node extreme vectors. We further prove and construct \\emph{sharp} all-$\\alpha$ corrections for the separate three-node chain, including exact globally feasible link parameters and an explicit polynomial positivity certificate.',1)
old=old.replace('The classical Kirchhoff transformation for temperature-dependent','The distinct three-node exact repair, certified below, supplements the earlier two-node construction without assuming that all fuzzy-marginal boxes are invalid. The classical Kirchhoff transformation for temperature-dependent',1)

section=r'''
\section{Attained three-node repair and whole-level feasibility}
\label{sec:3repair}
The three-node extrema in Section~4 quantify a conservation gap but do not on their own prove the closest jointly attainable state. Here $g_i=1$~W/K, so the temperature-sum deficit and energy gap have the same numerical value but different physical units. Let $\boldsymbol U(q)$ and $\boldsymbol L(q)$ be the \emph{exact} coordinate extrema derived from Eq.~\eqref{eq:3rational}, and define
\begin{equation}\label{eq:3gaps}
G_+(q)=\sum_{i=1}^3 U_i(q)-\tfrac{17}{5}(5+q),\quad
G_-(q)=\tfrac{17}{5}(5-q)-\sum_{i=1}^3 L_i(q).
\end{equation}
They are nonnegative and vanish at $q=0$. For any valid joint state $\boldsymbol t\in\mathcal S_\alpha$, we have $t_i\le U_i$ and $\sum_i t_i\le(17/5)(5+q)$, implying
\begin{equation}\label{eq:3lowerbounds}
\|\boldsymbol U-\boldsymbol t\|_1\ge G_+,
\qquad \|\boldsymbol U-\boldsymbol t\|_\infty\ge\tfrac{G_+}{3}.
\end{equation}
The lower extreme obeys the corresponding bounds $G_-$ and $G_-/3$. Because the ambient conductances are identical, choosing the \emph{source upper corner} yields $\sum_i t_i=\tfrac{17}{5}(5+q)$ for any allowed links. Thus the $L_1$ lower bound is automatically attained by any such state, and analogously for the lower source corner.

A nontrivial question is whether the three coordinate reductions can be made \emph{equal}, so that the $L_\infty$ bound is also sharp. We explicitly construct
\begin{equation}\label{eq:3balanced}
\boldsymbol t_+(q)=\boldsymbol U(q)-\tfrac{G_+(q)}{3}\boldsymbol1,\quad
\boldsymbol t_-(q)=\boldsymbol L(q)+\tfrac{G_-(q)}{3}\boldsymbol1.
\end{equation}
The accepted source vectors are $\boldsymbol Q_\pm=(5\pm q)(2,1,2/5)^\top$. For $s\in\{+,-\}$, define the link choices from the nodal flux balances, without any numerical optimizer:
\begin{equation}\label{eq:3linkstar}
a_s^*(q)=\frac{Q_{s,1}-t_{s,1}}{t_{s,1}-t_{s,2}},\qquad
b_s^*(q)=\frac{t_{s,3}-Q_{s,3}}{t_{s,2}-t_{s,3}}.
\end{equation}
Here the source notation $Q_{s,i}$ distinguishes the two corners from the earlier $Q_i$. Exact substitution proves $K(a_s^*,b_s^*)\boldsymbol t_s=\boldsymbol Q_s$. The remaining admission condition is $1-q/5\le a_s^*,b_s^*\le1+q/5$ at \emph{every} $q\in[0,1]$. After bringing all eight differences to rational form, each numerator is $qP(q)$ up to the constant sign of its denominator. All eight denominators have no zeros on $[0,1]$; each normalized polynomial satisfies the elementary whole-domain certificate
\begin{equation}\label{eq:3polygate}
P(q)\ge P(0)+\sum_{k\ge1}\min(0,P_k)>0,\quad
P(q)=\sum_kP_kq^k,\quad 0\le q\le1.
\end{equation}
The eight strictly positive integer lower bounds are $8{,}810{,}766$, $10{,}910{,}600$, $6{,}698{,}604$, $8{,}664{,}400$, $9{,}796{,}532$, $9{,}817{,}200$, $7{,}812{,}802$, and $7{,}423{,}800$ in the ordered machine-readable certificate. Consequently both balanced-repair states are admitted at all levels, establishing the \emph{attained} identities
\begin{equation}\label{eq:3sharp}
\boxed{\ d_{1,\pm}(q)=\frac{G_\pm(q)}{g_0^*},\quad d_{\infty,\pm}(q)=\frac{G_\pm(q)}{3g_0^*},\quad g_0^*=1\ \mathrm{W/K}.\ }
\end{equation}
Equivalently, with the numeric ambient conductance $1$~W/K, the values of $d_{1,\pm}$ and $d_{\infty,\pm}$ in K equal $G_\pm$ and $G_\pm/3$ in W numerically. This is a \emph{case-specific certificate}, not a universal theorem claiming equalization in every network.
\begin{table}[t]
\centering\small
\caption{Three-node support-level $q=1$ nearest jointly feasible corrections; units in K, and link conductances in W/K. All fractions are exact.}
\label{tab:3repair}
\begin{tabular}{@{}lcccc@{}}\toprule
Corner & $d_1$ (K) & $d_\infty$ (K) & $a^*$ (W/K) & $b^*$ (W/K) \\\midrule
Upper & $1855736/2541891$ & $1855736/7625673$ & $50692424/52306485$ & $39306238/40711005$ \\
Lower & $3695632/7625673$ & $3695632/22877019$ & $51261874/51768315$ & $39906968/40078995$ \\
\bottomrule\end{tabular}
\end{table}
\begin{figure}[t]
\centering\includegraphics[width=.80\linewidth]{Figure_R6_ThreeNode_SharpRepair.pdf}
\caption{Exact sharp $L_1$ and $L_\infty$ corrections for the original three-node synthetic network, derived analytically for every $\alpha$ and displayed as smooth curves. The curves are not solver-performance data.}\label{fig:3repair}
\end{figure}
A separate exact-rational replay verifies the two reconstructed states at 101 levels (202 local-balance substitutions). The continuum result rests on the eight polynomial inequalities and symbolic identities, \emph{not} on grid interpolation. The underlying projection-versus-joint-feasibility issue is established prior art; the new calculation provides a bounded auditable diagnostic for this published working benchmark.
'''
old=old.replace('\\section{Conditional temperature-dependent conductivity closure}', section+'\n\\section{Conditional temperature-dependent conductivity closure}',1)
old=old.replace('The matched-marginal two-node case now admits sharp exact $L_1$ and $L_\\infty$ joint-repair certificates for the entire membership interval.',
r'The matched-marginal two-node case and the independent-link three-node case now admit \emph{separate}, sharp, attained $L_1$ and $L_\\infty$ repair certificates for the entire membership interval. The three-node result uses an exact equal-component repair and eight whole-domain polynomial bounds.',1)
old=old.replace('The additional R5 NIST-Table-I1 audit script',
'The new R6 three-node symbolic verifier certifies eight all-level link-admissibility inequalities, with independent rational checks of 202 reconstructed thermal states. The additional R5 NIST-Table-I1 audit script',1)
# In body equations, use explicit ambient conductance g0*=1 W/K.

old=old.replace(chr(92)*2+'infty',chr(92)+'infty')
for fig in ('Figure_2node_Energy_Balance_Gap','Figure_NIST_MatchedMarginal_Dependence','Figure_NIST_2014_TableI1_SourceResiduals'):
    old=old.replace(fig+'.pdf',fig+'.png')
out=root/'ICMPC2027_R6_ThreeNode_SharpRepair_SourceAudited_WORKING.tex'
out.write_text(old,encoding='utf-8')
print('R6 document',out.name,'chars',len(old),'figure generated')
