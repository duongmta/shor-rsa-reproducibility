#!/usr/bin/env python3
"""Regenerate every number used in the paper's tables as CSV files in tables/.

    python paper_tables.py            # all parts whose inputs exist

  table2_per_base.csv / table2_grouped.csv   exact r, i, p(N,a), p_rand, G, m95   (Table 2, Fig. 2)
  table3_grouped.csv / table3_scaling.csv    logical gate counts before/after      (Table 3, Fig. 3)
  table4_compiled_all.csv                    Qiskit transpiled metrics             (Table 4)
  table5_aer_all.csv                         Qiskit Aer 1,024-shot statistics       (Table 5)
  independent_summary.csv                    independent statevector runs          (Section 5.4, Fig. 4b)
"""
import csv, glob, json, os, sys
from collections import defaultdict
from math import gcd, ceil, log
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import order, identity_blocks, success_probability, recovers_factor  # noqa: E402

OUT = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)


def write(name, rows, fields=None):
    if not rows:
        return
    fields = fields or list(rows[0].keys())
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote tables/{name} ({len(rows)} rows)")


def analytic():
    per, groups = [], []
    for N in (15, 21, 35):
        n = N.bit_length(); t = 2 * n; Q = 2 ** t
        g = defaultdict(list)
        for a in range(2, N):
            if gcd(a, N) != 1:
                continue
            r = order(a, N)
            p, _, _ = success_probability(N, a, t)
            prand = sum(recovers_factor(y, Q, a, N) for y in range(Q)) / Q
            rec = {"N": N, "a": a, "n": n, "t": t, "r": r, "i": identity_blocks(N, a, t),
                   "p": round(p, 6), "p_rand": round(prand, 6),
                   "G": round(p / prand, 3) if prand > 0 else "",
                   "m95": ceil(log(0.05) / log(1 - p)) if p > 0 else ""}
            per.append(rec)
            g[(r, round(p, 4), round(prand, 4))].append(rec)
        for k in sorted(g, key=lambda k: (k[0], -k[1])):
            L = g[k]
            groups.append({"N": N, "bases": " ".join(str(x["a"]) for x in L), "r": k[0], "i": L[0]["i"],
                           "p": f"{k[1]:.4f}", "p_rand": f"{k[2]:.4f}",
                           "G": f"{k[1] / k[2]:.1f}" if k[2] > 0 else "", "m95": L[0]["m95"]})
    write("table2_per_base.csv", per)
    write("table2_grouped.csv", groups)
    json.dump({"per_base": [{"N": r["N"], "a": r["a"], "r": r["r"], "id": r["i"], "p": r["p"],
                             "prand": r["p_rand"]} for r in per]},
              open(os.path.join(HERE, "paper_analytic.json"), "w"), indent=1)
    for N in (15, 21, 35):
        print(f"  mean p(N,a) over bases, N={N}: {np.mean([r['p'] for r in per if r['N'] == N]):.4f}")


def table3():
    path = os.path.join(HERE, "logical_counts.csv")
    if not os.path.exists(path):
        print("  skip Table 3: run local_tables.py first"); return
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    g = defaultdict(list)
    for r in rows:
        g[(int(r["N"]), int(r["order"]), int(r["identity_blocks"]))].append(r)
    out = []
    for (N, r_, i_), L in sorted(g.items()):
        rec = {"N": N, "r": r_, "i": i_, "bases": len(L) // 3}
        for fam, tag in (("beauregard", "A"), ("cuccaro", "B"), ("montgomery", "C")):
            rr = [x for x in L if x["family"] == fam]
            b = np.mean([int(x["before_gates"]) for x in rr]); a = np.mean([int(x["after_gates"]) for x in rr])
            bs = [int(x["before_gates"]) for x in rr]
            rec[f"{tag}_before"] = round(b); rec[f"{tag}_after"] = round(a)
            rec[f"{tag}_reduction_%"] = round(100 * (1 - a / b), 1)
            rec[f"{tag}_spread_%"] = round(100 * (max(bs) - min(bs)) / b, 2)
        out.append(rec)
    write("table3_grouped.csv", out)
    sc = []
    for fam, tag in (("beauregard", "A"), ("cuccaro", "B"), ("montgomery", "C")):
        means = [np.mean([int(x["before_gates"]) for x in rows if x["family"] == fam and int(x["N"]) == N])
                 for N in (15, 21, 35)]
        slope = np.polyfit(np.log([4, 5, 6]), np.log(means), 1)[0]
        sc.append({"family": tag, "mean_before_N15": round(means[0]), "mean_before_N21": round(means[1]),
                   "mean_before_N35": round(means[2]), "fitted_exponent": round(slope, 2)})
    write("table3_scaling.csv", sc)
    red = [float(x["gate_reduction_%"]) for x in rows if int(x["identity_blocks"]) > 0]
    print(f"  reduction range over bases with identity blocks: {min(red)}-{max(red)} %")


def table4():
    path = os.path.join(HERE, "results", "counts.jsonl")
    if not os.path.exists(path):
        print("  skip Table 4: run run_qiskit.py counts first"); return
    recs = [json.loads(l) for l in open(path, encoding="utf-8")]
    before = {(r["family"], r["N"], r["a"]): r for r in recs if r["stage"] == "before"}
    out = []
    for r in sorted(recs, key=lambda r: (r["family"], r["N"], r["a"], r["stage"] != "before")):
        L, T = r["logical"], r["transpiled"]
        b = before.get((r["family"], r["N"], r["a"]))
        out.append({"family": r["family"], "N": r["N"], "a": r["a"], "stage": r["stage"], "r": r["order"],
                    "i": r["identity_blocks"], "qubits": r["qubits"], "logical_ops": L["size"],
                    "logical_depth": L["depth"], "basis_ops": T["size"], "cx": T["cx"], "depth": T["depth"],
                    "rz_clifford": T["rz_clifford"], "rz_t_like": T["rz_t_like"],
                    "rz_arbitrary": T["rz_arbitrary"],
                    "cx_reduction_%": round(100 * (1 - T["cx"] / b["transpiled"]["cx"]), 1) if b else "",
                    "depth_reduction_%": round(100 * (1 - T["depth"] / b["transpiled"]["depth"]), 1) if b else "",
                    "transpile_s": r["transpile_s"]})
    write("table4_compiled_all.csv", out)


def table5():
    path = os.path.join(HERE, "results", "simulate.jsonl")
    if not os.path.exists(path):
        print("  skip Table 5: run run_qiskit.py simulate first"); return
    recs = [json.loads(l) for l in open(path, encoding="utf-8")]
    out = []
    for r in sorted(recs, key=lambda r: (r["family"], r["N"], r["a"], r["stage"] != "before")):
        top = sorted(r["counts"].items(), key=lambda kv: -kv[1])[:6]
        out.append({"family": r["family"], "N": r["N"], "a": r["a"], "stage": r["stage"], "r": r["order"],
                    "qubits": r["qubits"], "method": r["method"], "seed": r["seed_simulator"],
                    "shots": r["shots"], "p_exact": round(r["p_exact"], 4), "success_count": r["success_count"],
                    "p_hat": round(r["p_hat"], 4), "wilson_lo": round(r["wilson95"][0], 4),
                    "wilson_hi": round(r["wilson95"][1], 4), "p_exact_in_ci": r["p_exact_in_ci"],
                    "tvd": round(r["tvd"], 4), "chi2": round(r["chi2"], 2), "dof": r["dof"],
                    "chi2_p": round(r["chi2_p"], 3), "off_support": r["off_support"], "wall_s": r["wall_s"],
                    "top_outcomes(y:count)": " ".join(f"{y}:{c}" for y, c in top)})
    write("table5_aer_all.csv", out)
    print(f"  Aer runs: {len(out)}, all p_exact in CI: {all(r['p_exact_in_ci'] for r in out)}, "
          f"total off-support shots: {sum(r['off_support'] for r in out)}")


def independent():
    files = sorted(glob.glob(os.path.join(HERE, "results_independent", "*.jsonl")))
    if not files:
        print("  skip independent runs: none found"); return
    out = []
    for fpath in files:
        for l in open(fpath, encoding="utf-8"):
            r = json.loads(l)
            out.append({"source_file": os.path.basename(fpath), "family": r["family"], "N": r["N"], "a": r["a"],
                        "qubits": r["qubits"], "gates": r["gates"], "r": r["order"], "i": r["identity_blocks"],
                        "max_dev": f"{r['max_dev']:.2e}", "p_exact": round(r["p_exact"], 4),
                        "success_count": r["success_count"], "p_hat": round(r["p_hat"], 4),
                        "wilson_lo": round(r["wilson95"][0], 4), "wilson_hi": round(r["wilson95"][1], 4),
                        "in_ci": r["wilson95"][0] <= r["p_exact"] <= r["wilson95"][1],
                        "tvd": round(r["tvd"], 4), "chi2_p": round(r["chi2_p"], 3),
                        "off_support": r["off_support"], "wall_s": r["wall_s"]})
    write("independent_summary.csv", out)
    for src in sorted({o["source_file"] for o in out}):
        rr = [o for o in out if o["source_file"] == src]
        print(f"  {src}: {len(rr)} circuits, max_dev {max(float(o['max_dev']) for o in rr):.1e}, "
              f"in CI {sum(o['in_ci'] for o in rr)}/{len(rr)}, off-support {sum(o['off_support'] for o in rr)}")


if __name__ == "__main__":
    for step in (analytic, table3, table4, table5, independent):
        print(f"[{step.__name__}]")
        step()
