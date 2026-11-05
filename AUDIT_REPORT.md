# Audit of the reproducibility package

Audit date: **5 November 2026**. Result: **PASS for the checks listed below**.

This package was assembled from `shor_reversible/` (original generators and evidence) and the two scientific programs plus eight data/environment files in `revision_20261104/`. All 18 original scientific Python files and all other selected original evidence files are byte-identical to their source. The 64 copied files were also checked against the prior reproducibility ZIP, whose file manifest passed. The current manuscript is identified by its hash in `MANUSCRIPT_CONTEXT.json`; no earlier manuscript was used as the numerical target.

## Checks performed afresh

- Rebuilt and recompiled **all 246 circuits**, using Qiskit 2.2.3. Every recorded logical and compiled metric matched, including gate counts, depths and rotation classes.
- Recomputed exact finite-QPE probabilities and uniform-reference acceptance for all **41 pairs**. Recovered the order-two acceptance numerators 123, 523 and 2337.
- Recalculated acceptance, Wilson intervals, TVD and chi-square statistics from **all 144 archived Aer histograms**. Checked current Tables 1–6 against their data sources.
- Regenerated all **seven historical CSV summary tables**; their rows matched the archive.
- Ran the simulator sanity check and verified **132 controlled-multiplier blocks**: 44 constants for each family, with both control values and every x<N. The checks establish their stated basis-input mapping and cleanup, not full phase-coherent equivalence.
- Repeated all small-instance post-processing/precision cases and the nine larger analytic-distribution sampling cases. Both sensitivity JSON files and the spread calculation matched exactly.
- Repeated the four A-family sampling diagnostics. Before/after counts matched their archived values with truncation enabled, and the archived equal pair was recovered with truncation disabled. Execution metadata/timing were not required to be identical.
- Repeated all 12 compiler comparisons, B's two arithmetic Clifford+T decompositions, and all three terminal inverse-QFT settings. Their reported numerical results matched exactly.

The detailed machine-readable reports and fresh logs are in `validation/`. Individual simulator/compiler runtimes are machine-dependent; historic times in the paper were preserved.

## Scope of this audit

The full 144-run Aer sweep and the full 25-configuration independent statevector sweep were **not rerun** during this packaging audit. Their archived outputs, generator provenance and coverage were checked. New Aer calculations were confined to the sampling diagnostic and the approximate-QFT tests. The repository includes the programs and fresh-output wrapper needed to repeat the longer sweeps separately.

There are 246 resource configurations, not 246 archived Aer simulations. The five independent JSONL files include cloud/local duplicates; they represent 25 distinct after configurations. One separate before check is preserved as a text log. The 10–12-bit extension is analytic-distribution sampling, not full-circuit simulation.

## Files selected and excluded

Included: unchanged scientific source, raw JSONL/JSON/CSV data, relevant historical logs, original plotting programs, the current manuscript's image assets, pinned numerical dependencies, evidence maps, fresh-run wrappers and verification scripts.

Excluded: Word/PDF manuscripts, reviewer comments and responses, obsolete supplementary-appendix documents, document-formatting programs, old instructions with inaccurate version/coverage statements, the old batch launcher, caches, duplicate extracted packages, and unused auxiliary plot exports. The archived original inputs remain in their original locations; this release was made by copying, not moving or deleting them.

Targeted scans found no matches for private-key, GitHub-token, AWS-key or common credential-assignment patterns in the distributed text files. This is a limited automated scan, not an absolute security guarantee.

Git is initialized locally solely to simplify use with GitHub Desktop. No remote repository, public URL or DOI was created. Publication and the final commit identity remain with the authors. Instructions are in `HUONG_DAN_GITHUB_VI.md`.
