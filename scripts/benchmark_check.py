# scripts/benchmark_check.py
import pandas as pd
from mpx_dpcr.evaluate import Oligo, evaluate_set, pass_fail_summary

# === Benchmark oligos from Hindson et al. (Anal Chem 2011) ===
# NOTE: For probes, include the plain DNA sequence (no fluor/quencher text).
# MRGPRX1 (FAM-labeled probe in paper)
MRGPRX1_F = "TTAAGCTTCATCAGTATCCCCCA"
MRGPRX1_R = "CAAAGTAGGAAAACATCATCACAGGA"
MRGPRX1_PROBE = "ACCATCTCTAAAATCCT"  # reported as 6FAM-...-MGBNFQ

# RPP30 (VIC-labeled probe in paper)
RPP30_F = "GATTTGGACCTGCGAGCG"
RPP30_R = "GCGGCTGTCTCCACAAGT"
RPP30_PROBE = "CTGACCTGAAGGCTCT"  # reported as VIC-...-MGBNFQ

oligos = [
    Oligo("MRGPRX1_F", MRGPRX1_F, "primer_f"),
    Oligo("MRGPRX1_R", MRGPRX1_R, "primer_r"),
    Oligo("MRGPRX1_probe", MRGPRX1_PROBE, "probe"),
    Oligo("RPP30_F", RPP30_F, "primer_f"),
    Oligo("RPP30_R", RPP30_R, "primer_r"),
    Oligo("RPP30_probe", RPP30_PROBE, "probe"),
]

per_oligo, cross, violations = evaluate_set(oligos)
summary = pass_fail_summary(violations)

per_oligo.to_csv("benchmark_per_oligo.csv", index=False)
cross.to_csv("benchmark_cross_matrix.csv")
violations.to_csv("benchmark_violations.csv", index=False)

print("PASS:", summary["passes_cutoff"], "Violations:", summary["n_violations"])
print("Wrote: benchmark_per_oligo.csv, benchmark_cross_matrix.csv, benchmark_violations.csv")
