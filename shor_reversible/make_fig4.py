"""Figure 4 -> figures/Hinh_4.png: (a) Aer N=15,a=7 (results/simulate.jsonl), (b) independent runs."""
import json, os, glob, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6, "axes.axisbelow": True})
B, O, G, GR = "#2a78d6", "#eb6834", "#1baf7a", "#9a9a9a"
os.makedirs("figures", exist_ok=True)
aer = [json.loads(l) for l in open("results/simulate.jsonl", encoding="utf-8")]
aer = [q for q in aer if q["N"] == 15 and q["a"] == 7]
def _load(names):
    for n in names:
        p = os.path.join("results_independent", n)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return [json.loads(l) for l in open(p, encoding="utf-8")]
    return []
# prefer runs made on this computer (local_*.jsonl); fall back to the cloud runs used in the paper
cloudA = _load(["independent_A.jsonl"])
A15 = _load(["local_A15.jsonl"]) or [r for r in cloudA if r["N"] == 15]
A21 = _load(["local_A21.jsonl"]) or [r for r in cloudA if r["N"] == 21]
B15 = _load(["local_B15.jsonl", "independent_B15.jsonl"])
sw = A15 + A21 + B15
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 3.0), gridspec_kw={"width_ratios": [1, 2.1]})
ys = [0, 64, 128, 192]; x = np.arange(4); w = 0.2
ax1.bar(x - 1.5 * w, [0.25] * 4, w, color=GR, label="Exact")
for k, (fam, col, lab) in enumerate((("beauregard", B, "A"), ("cuccaro", O, "B"), ("montgomery", G, "C"))):
    r = [q for q in aer if q["family"] == fam and q["stage"] == "after"][0]
    f = [r["counts"].get(str(y), 0) / r["shots"] for y in ys]
    ax1.bar(x + (k - 0.5) * w, f, w, color=col, label=f"Aer, family {lab}")
ax1.set_xticks(x); ax1.set_xticklabels([str(y) for y in ys]); ax1.set_xlabel("Outcome y"); ax1.set_ylabel("Probability / frequency")
ax1.set_ylim(0, 0.45); ax1.set_yticks([0, 0.1, 0.2, 0.3]); ax1.legend(frameon=False, fontsize=6.3, ncol=2, loc="upper right", bbox_to_anchor=(1.02, 1.02), columnspacing=0.8, handlelength=1.2)
ax1.set_title("(a) Qiskit Aer, N = 15, a = 7", fontsize=8.5, loc="left")
groups = [("beauregard", 15, B, "o", "A, N = 15"), ("cuccaro", 15, O, "s", "B, N = 15"), ("beauregard", 21, B, "D", "A, N = 21")]
pos = 0; ticks = []; labs = []
for fam, N, col, mk, lab in groups:
    rr = sorted([r for r in sw if r["family"] == fam and r["N"] == N], key=lambda r: r["a"])
    xs = np.arange(pos, pos + len(rr))
    ax2.hlines([r["p_exact"] for r in rr], xs - 0.38, xs + 0.38, color="#222222", lw=1.1, zorder=2)
    ph = np.array([r["p_hat"] for r in rr]); lo = ph - np.array([r["wilson95"][0] for r in rr]); hi = np.array([r["wilson95"][1] for r in rr]) - ph
    ax2.errorbar(xs, ph, yerr=[lo, hi], fmt=mk, ms=3.8, color=col, lw=0.9, capsize=0, label=lab, zorder=3,
                 mfc=col if fam == "beauregard" else "white")
    ticks += list(xs); labs += [str(r["a"]) for r in rr]
    pos += len(rr) + 1
ax2.set_xticks(ticks); ax2.set_xticklabels(labs, fontsize=6.3, rotation=90); ax2.set_xlabel("Base a")
ax2.set_ylabel("Factor-recovery fraction"); ax2.set_ylim(-0.03, 0.74); ax2.set_yticks(np.arange(0, 0.61, 0.1))
ax2.plot([], [], color="#222222", lw=1.1, label="Exact p(N, a)")
ax2.legend(frameon=False, fontsize=6.5, ncol=4, loc="upper center")
ax2.set_title("(b) Independent statevector, 1,024 samples per base", fontsize=8.5, loc="left")
fig.tight_layout(w_pad=1.0); fig.savefig("figures/Hinh_4.png", dpi=300)
print("fig4 done", len(sw))
