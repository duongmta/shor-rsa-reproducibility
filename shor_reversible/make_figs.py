"""Figures 1-3 of the paper -> figures/Hinh_1..3.png (needs paper_analytic.json from paper_tables.py
and logical_counts.csv from local_tables.py)."""
import json, csv, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6,
                     "mathtext.fontset": "dejavuserif"})
B, O, G = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#222222", "#6b6b6b"
os.makedirs("figures", exist_ok=True)
A = json.load(open("paper_analytic.json"))

# ---------------- Fig. 1: architecture schematic ----------------
fig, ax = plt.subplots(figsize=(6.3, 3.3)); ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")
ys = [52, 47.5, 43, 36]            # phase lines P0,P1,P2,...,P_{t-1}
labels = [r"$P_0$", r"$P_1$", r"$P_2$", r"$P_{t-1}$"]
wy = 27                            # work register
for y, l in zip(ys, labels):
    ax.plot([12, 88], [y, y], color=INK, lw=0.8, zorder=1); ax.text(0.5, y, l + r"$|0\rangle$", va="center", fontsize=8)
ax.text(15.2, 39.8, r"$\vdots$", fontsize=9, ha="center")
ax.add_patch(Rectangle((12, wy - 3.2), 76, 6.4, fc="#f3f3f3", ec="none", zorder=0))
ax.plot([12, 88], [wy, wy], color=INK, lw=2.2, zorder=1)
ax.text(0.5, wy, r"$|1\rangle_n|0\cdots\rangle$", va="center", fontsize=8)
ax.text(47, wy - 5.2, "work register + arithmetic ancillas (family-dependent)", ha="center", fontsize=7, color=MUTED)
for y in ys:
    ax.add_patch(Rectangle((13.5, y - 1.7), 3.4, 3.4, fc="white", ec=INK, lw=0.8, zorder=3)); ax.text(15.2, y, "H", zorder=4, ha="center", va="center", fontsize=7)
blocks = [(19, 0, r"$U_{c_0}$", False), (31, 1, r"$U_{c_1}$", False), (43, 2, r"$U_{c_2}$", True), (60, 3, r"$U_{c_{t-1}}$", True)]
for x, j, lab, ident in blocks:
    ax.plot([x + 4, x + 4], [ys[j], wy + 3.2], color=MUTED if ident else INK, lw=0.8, ls=(0, (2, 1.5)) if ident else "-")
    ax.plot(x + 4, ys[j], "o", ms=3.5, color=MUTED if ident else INK, zorder=4)
    ax.add_patch(Rectangle((x, wy - 3.2), 8, 6.4, fc="white" if not ident else "#e9e9e9", ec=MUTED if ident else INK,
                           lw=0.8, hatch="////" if ident else None, zorder=3))
    ax.text(x + 4, wy, lab, ha="center", va="center", fontsize=8, color=MUTED if ident else INK, zorder=4,
            bbox=dict(fc="white", ec="none", pad=0.6) if ident else None)
ax.text(53, wy, r"$\cdots$", ha="center", va="center"); ax.text(53, 44, r"$\cdots$", ha="center")
ax.add_patch(Rectangle((73, 34), 7, 20, fc="white", ec=INK, lw=0.8, zorder=3)); ax.text(76.5, 44, r"QFT$^{\dagger}$", ha="center", va="center", rotation=90, fontsize=8, zorder=4)
for y in ys:
    ax.add_patch(Rectangle((83, y - 1.7), 4, 3.4, fc="white", ec=INK, lw=0.8, zorder=3)); ax.text(85, y, "M", zorder=4, ha="center", va="center", fontsize=6.5)
ax.text(55, 17.8, r"$c_j=a^{2^j}\ \mathrm{mod}\ N$;  hatched blocks have $c_j=1$ (removed in the after-circuit)", ha="center", fontsize=7.2, color=MUTED)
fams = [(1, B, "A  Beauregard (Fourier)", "CMULT$(k)$ → CSWAP → CMULT$(k^{-1})^{\\dagger}$\nDraper $\\Phi$ADD + modular adder\n$2n+2$ work qubits"),
        (34.5, O, "B  Toffoli ripple-carry", "Cuccaro MAJ/UMA adder, VBE-style\nreduction, shared constant register\n$3n+3$ work qubits"),
        (68, G, "C  Montgomery", "REDC: $+k'x_i$, parity $\\to q_i$, $+q_iN$, shift;\ncopy $y$, Bennett uncompute\n$4n+3$ work qubits")]
for x, col, title, body in fams:
    ax.add_patch(FancyBboxPatch((x, 0.5), 31, 14.5, boxstyle="round,pad=0.2,rounding_size=1.2", fc="white", ec=col, lw=1.3))
    ax.text(x + 1.5, 12.6, title, fontsize=7.6, color=INK, weight="bold", va="center")
    ax.text(x + 1.5, 6.2, body, fontsize=6.6, color=INK, va="center", linespacing=1.35)
fig.tight_layout(pad=0.2); fig.savefig("figures/Hinh_1.png", dpi=300); plt.close(fig)

# ---------------- Fig. 2: base dependence (exact) ----------------
fig, axes = plt.subplots(3, 1, figsize=(6.3, 5.6))
for ax, N in zip(axes, (15, 21, 35)):
    rows = [r for r in A["per_base"] if r["N"] == N]
    x = np.arange(len(rows))
    ax.bar(x, [r["p"] for r in rows], 0.62, color=B, label="Exact finite-QPE $p(N,a)$", zorder=2)
    ax.plot(x, [r["prand"] for r in rows], "D", ms=4.2, color=O, mec="white", mew=0.6, label="Uniform-random phase, same rule", zorder=3)
    ax.set_xticks(x); ax.set_xticklabels([f"{r['a']}\n{r['r']}\n{r['id']}" for r in rows], fontsize=6.3)
    ax.set_ylim(0, 0.66); ax.set_ylabel("Probability"); ax.set_xlim(-0.7, len(rows) - 0.3)
    ax.text(0.005, 0.97, f"N = {N}  (t = {2*N.bit_length()})", transform=ax.transAxes, va="top", fontsize=8, weight="bold")
    ax.text(-0.075, -0.2, "a\nr\nid", transform=ax.transAxes, fontsize=6.3, va="center", ha="right", color=MUTED)
axes[0].legend(frameon=False, fontsize=7.2, loc="upper right", ncol=2)
fig.tight_layout(h_pad=0.6); fig.savefig("figures/Hinh_2.png", dpi=300); plt.close(fig)

# ---------------- Fig. 3: gate reduction vs kept-block fraction ----------------
rows = list(csv.DictReader(open("logical_counts.csv", encoding="utf-8")))
fig, ax = plt.subplots(figsize=(3.9, 2.9))
ax.plot([0, 1], [0, 1], color="#bbbbbb", lw=1, zorder=1, label="Ideal: ratio = kept fraction")
for fam, col, mk, lab in (("beauregard", B, "o", "A  Beauregard"), ("cuccaro", O, "s", "B  Toffoli"), ("montgomery", G, "^", "C  Montgomery")):
    rr = [r for r in rows if r["family"] == fam]
    xs = [(2 * int(r["N"]).bit_length() - int(r["identity_blocks"])) / (2 * int(r["N"]).bit_length()) for r in rr]
    ys_ = [int(r["after_gates"]) / int(r["before_gates"]) for r in rr]
    ax.scatter(xs, ys_, s=22, marker=mk, facecolor="none" if fam != "beauregard" else col, edgecolor=col, lw=1, label=lab, zorder=3)
ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
ax.set_xlabel("Fraction of non-identity blocks, $(t-i)/t$"); ax.set_ylabel("Logical gates after / before")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig("figures/Hinh_3.png", dpi=300); plt.close(fig)
print("figs 1-3 done")
