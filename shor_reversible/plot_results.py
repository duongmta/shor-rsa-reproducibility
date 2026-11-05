#!/usr/bin/env python3
"""Figures from results/*.jsonl (needs matplotlib).

  python plot_results.py            -> results/fig_*.png (300 dpi)
Fig. hist_<fam>_<N>_<a>: exact QPE distribution vs Aer frequencies (1,024 shots)
Fig. success_<N>:        p_hat (Wilson 95% CI) vs exact p(N,a) for every base a, 3 families
Fig. counts_<N>:         basis-level CX count before/after identity-block removal, per base
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import order, qpe_distribution   # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                 # noqa: E402

OUT = os.path.join(HERE, "results")
COL = {"beauregard": "#2a78d6", "cuccaro": "#eb6834", "montgomery": "#1baf7a"}
LAB = {"beauregard": "Beauregard (QFT adder)", "cuccaro": "Toffoli (Cuccaro adder)",
       "montgomery": "Montgomery"}
MARK = {"beauregard": "o", "cuccaro": "s", "montgomery": "^"}
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                     "grid.linewidth": 0.6})


def load(name):
    p = os.path.join(OUT, name)
    return [json.loads(l) for l in open(p, encoding="utf-8")] if os.path.exists(p) else []


def fig_hist(rec):
    N, a, t = rec["N"], rec["a"], 2 * rec["N"].bit_length()
    Q = 2 ** t
    P = qpe_distribution(order(a, N), Q)
    obs = np.zeros(Q)
    for y, c in rec["counts"].items():
        obs[int(y)] += c
    f = obs / rec["shots"]
    fig, ax = plt.subplots(figsize=(6.3, 2.6))
    y = np.arange(Q)
    ax.vlines(y - Q * 0.002, 0, P, color="#2a78d6", lw=2, label="Exact finite-QPE probability")
    ax.vlines(y + Q * 0.002, 0, f, color="#eb6834", lw=2, label=f"Qiskit Aer frequency ({rec['shots']} shots)")
    ax.set_xlim(-Q * 0.02, Q * 1.02)
    ax.set_xlabel("Measured phase integer y"); ax.set_ylabel("Probability")
    ax.set_title(f"{LAB[rec['family']]}, N = {N}, a = {a}, t = {t}   "
                 f"(TVD = {rec['tvd']:.4f}, off-support shots = {rec['off_support']})", fontsize=8.5)
    ax.set_ylim(0, max(P.max(), f.max()) * 1.4)
    ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_hist_{rec['family']}_{N}_{a}.png"), dpi=300)
    plt.close(fig)


def fig_success(sim, N):
    rows = [r for r in sim if r["N"] == N and r["stage"] == "after"]
    if not rows: return
    bases = sorted({r["a"] for r in rows})
    x = {a: i for i, a in enumerate(bases)}
    fig, ax = plt.subplots(figsize=(6.3, 2.8))
    exact = {r["a"]: r["p_exact"] for r in rows}
    ax.hlines([exact[a] for a in bases], [x[a] - 0.4 for a in bases], [x[a] + 0.4 for a in bases],
              color="#333333", lw=1.2, label="Exact p(N, a)")
    fams = [f for f in COL if any(r["family"] == f for r in rows)]
    for k, fam in enumerate(fams):
        rr = [r for r in rows if r["family"] == fam]
        off = (k - (len(fams) - 1) / 2) * 0.22
        xs = [x[r["a"]] + off for r in rr]
        ph = [r["p_hat"] for r in rr]
        lo = [r["p_hat"] - r["wilson95"][0] for r in rr]
        hi = [r["wilson95"][1] - r["p_hat"] for r in rr]
        ax.errorbar(xs, ph, yerr=[lo, hi], fmt=MARK[fam], ms=4, color=COL[fam], lw=1,
                    capsize=0, label=LAB[fam])
    ax.set_xticks(range(len(bases)))
    ax.set_xticklabels([f"{a}\n(r={order(a, N)})" for a in bases], fontsize=7)
    ax.set_xlabel("Base a (order r)"); ax.set_ylabel("Factor-producing fraction")
    ax.set_ylim(-0.03, 0.62)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_success_{N}.png"), dpi=300)
    plt.close(fig)


def fig_counts(cnt, N):
    rows = [r for r in cnt if r["N"] == N]
    if not rows: return
    fams = [f for f in COL if any(r["family"] == f for r in rows)]
    bases = sorted({r["a"] for r in rows})
    fig, axes = plt.subplots(1, len(fams), figsize=(6.3, 2.5), sharey=True, squeeze=False)
    for ax, fam in zip(axes[0], fams):
        get = {(r["a"], r["stage"]): r["transpiled"]["cx"] for r in rows if r["family"] == fam}
        xi = np.arange(len(bases))
        b = [get.get((a, "before"), np.nan) for a in bases]
        a_ = [get.get((a, "after"), np.nan) for a in bases]
        ax.bar(xi - 0.2, b, 0.38, color="#b7d3f6", label="All t blocks")
        ax.bar(xi + 0.2, a_, 0.38, color=COL[fam], label="Identity blocks removed")
        ax.set_title(LAB[fam], fontsize=8)
        ax.set_xticks(xi); ax.set_xticklabels(bases, fontsize=6.5); ax.set_xlabel("Base a")
    axes[0][0].set_ylabel("CX count after transpilation")
    axes[0][0].legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_counts_{N}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    sim = load("simulate.jsonl"); cnt = load("counts.jsonl")
    for rec in sim:
        if rec["stage"] == "after":
            fig_hist(rec)
    for N in sorted({r["N"] for r in sim}):
        fig_success(sim, N)
    for N in sorted({r["N"] for r in cnt}):
        fig_counts(cnt, N)
    print("figures written to", OUT)
