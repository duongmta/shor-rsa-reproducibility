# Base-resolved Shor order-finding circuits

Code and recorded data supporting **Base-resolved evaluation of generic reversible circuits for Shor order finding on small RSA moduli**. This package is tied to the revised manuscript snapshot dated 5 November 2026 by the SHA-256 in `MANUSCRIPT_CONTEXT.json`. It is prepared locally for publication; no repository URL or DOI has been assigned by this package.

The experiment compares three reversible modular-multiplication implementations: A (`beauregard`, Fourier arithmetic), B (`cuccaro`, ripple-carry arithmetic), and C (`montgomery`, Montgomery reduction). The `before` circuit retains all controlled multipliers; `after` omits precisely those whose known multiplier constant is 1 modulo N. The circuit generator does not take the factors or multiplicative order as inputs.

## Scope and contents

| Evidence | Coverage | Location |
|---|---|---|
| Compiled resources | 41 (N,a) pairs × 3 families × 2 stages = **246 circuits** | `shor_reversible/results/counts.jsonl` |
| Logical resources | 123 rows, each containing before and after | `shor_reversible/logical_counts.csv` |
| Aer samples | **144 runs**: all 123 after circuits and 21 before circuits at N=15; 1,024 shots each | `shor_reversible/results/simulate.jsonl` |
| Independent statevector data | 25 distinct after circuits, with a separate before check in the logs | `shor_reversible/results_independent/`, `shor_reversible/logs/05_e2e_A15_a7_before.txt` |
| Alternative post-processing | 41 bases × 3 phase widths × 2 rules = 246 analytic records, not 246 additional circuits | `revision_20261104/postprocessing_small.json` |
| Larger moduli | Nine selected (N,a) pairs, 20,000 samples per distribution, two rules on the same samples | `revision_20261104/postprocessing_large.json` |
| Compiler, T resources and approximate QFT | Revision tests described in Section 4.3 | `revision_20261104/*` |

`BASES_41.csv` lists all pairs: 7 bases for N=15, 11 for N=21 and 23 for N=35. The 246 resource records do **not** mean that 246 full circuits were simulated by Aer.

`SOURCE_PROVENANCE.json` records the original relative location and hash of every copied computational file. All scientific Python files were copied without changes. `EVIDENCE_MAP.csv` maps the current paper's six tables and additional results to code and data. `MANIFEST_SHA256.json` detects changes to distributed files.

## Environment and installation

Recorded compiled/Aer and revision experiments used **Python 3.12.10, Qiskit 2.2.3 and Qiskit Aer 0.17.2** on Windows 11. The current reproduction environment also records NumPy 1.26.4, SciPy 1.15.3 and Matplotlib 3.10.7. Package versions are pinned in `requirements.txt`; the earlier independent cloud environment was not fully pinned.

In PowerShell, inside this repository:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -B check_release.py
```

The default check verifies hashes and coverage, rebuilds all **246 logical circuits**, recomputes exact probabilities for all 41 bases, recalculates statistics from all 144 archived histograms, and checks Tables 1–6 against the manuscript snapshot. It does not pretend to rerun archived quantum simulations.

To additionally **recompile all 246 circuits** and compare every stored logical and compiled metric:

```powershell
.\.venv\Scripts\python.exe -B check_release.py --compiled
```

Different compiler versions can change resource counts. A failing comparison should be investigated, rather than replacing the archived measurements with new ones.

## Fresh reproductions without overwriting evidence

The original programs have fixed output paths and some resume existing results. Use the wrapper below, which creates a fresh working copy under `runs/` and refuses to reuse an existing output directory:

```powershell
.\.venv\Scripts\python.exe -B reproduce.py tables
.\.venv\Scripts\python.exe -B reproduce.py analytic
.\.venv\Scripts\python.exe -B reproduce.py quantum
.\.venv\Scripts\python.exe -B reproduce.py aqft
```

`tables` regenerates the original seven CSV summaries from archived inputs and checks equality. `analytic` recomputes the small-base sensitivity analysis, Table 6 and within-group spreads. `quantum` repeats the seeded-sampling diagnostic, optimization levels 1/3 and the arithmetic-only B T metrics. `aqft` repeats the terminal inverse-QFT comparison at N=21,a=2. New JSON outputs can be compared with the corresponding archived files. Sampling and floating-point results may depend on the environment.

Additional modes are `counts` (fresh compilation), `functional` (simulator sanity and all 132 multiplier blocks), `aer` (all 144 Aer configurations), and `independent` (the 25 full-statevector after checks and one before check). Full simulation modes can take substantially longer and need more memory. For a second run, use a different directory, for example `--out runs/analytic_second`.

The `aer` mode first runs the six Table-5 configurations using their original automatic method selection, then fills the remaining configurations with MPS, matching the archived method choices. Do not substitute a uniform MPS rerun for the six original Table-5 runs.

## Interpretation and figure provenance

The basis-input checks require the expected state mapping and cleanup probability at least 1−10⁻⁹; they do not alone prove coherent operator equivalence. Independent deviations in the paper concern marginal probabilities, not complex amplitudes. Archived duplicate cloud/local runs are not additional distinct circuits. MPS is floating-point simulation, not exact arithmetic.

The 10–12-bit extension samples analytic QPE distributions; it does not build or simulate those complete quantum circuits. T-depth applies to one fixed decomposition of **B arithmetic only**. The A/C rotation-cost estimates are leading-order synthesis estimates, not completed fault-tolerant compilations.

`figures_manuscript/` contains the exact five image assets embedded in the current manuscript: the two parts of Fig. 1 and Figs. 2–4. They are mapped and hashed in `MANUSCRIPT_CONTEXT.json`. `shor_reversible/figures/` and the plotting scripts retain historical plot outputs; they need not match Word-resized/reformatted image bytes. The Fig. 1 internal-flow panel is an authored schematic, not experimental data.

The old project README, old batch launcher, Word-editing programs, reviewer documents, obsolete appendix document, caches and duplicate extracted packages are intentionally absent. Their obsolete version and coverage claims are superseded by this README. See `AUDIT_REPORT.md` for the checks actually performed on this release.
