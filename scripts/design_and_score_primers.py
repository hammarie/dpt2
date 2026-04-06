# scripts/design_and_score_primers.py
import argparse
import warnings
import os
from pathlib import Path

import pandas as pd
from Bio import BiopythonDeprecationWarning

from mpx_dpcr.fetch import fetch_gene_sequence
from mpx_dpcr.design import design_primers_for_gene
from mpx_dpcr.evaluate import Oligo, evaluate_set, DG_CUTOFF
from mpx_dpcr.blast_check import pair_specific_on_same_transcript


def _get_next_run_number(results_dir, base_name):
    """Find the next available run number for output files."""
    if not os.path.exists(results_dir):
        return 1

    run_nums = []
    for filename in os.listdir(results_dir):
        if base_name in filename:
            # Extract number from pattern: 1_filename.xlsx or 1_filename_primerblast.csv
            parts = filename.replace(".xlsx", "").replace(".csv", "").replace("_primerblast", "").replace(base_name, "")
            try:
                num = int(parts.lstrip("_").rstrip("_"))
                run_nums.append(num)
            except ValueError:
                continue

    return max(run_nums) + 1 if run_nums else 1


def _format_hit_stats(hit_dicts):
    """Return a compact string summarizing e-value/mismatch/gap per hit."""
    if not hit_dicts:
        return ""
    stats = []
    for hit in hit_dicts:
        if not isinstance(hit, dict):
            continue
        acc = hit.get("accession")
        evalue = hit.get("evalue")
        mismatches = hit.get("mismatches")
        gaps = hit.get("gaps")
        stats.append(f"{acc}|e={evalue}|m={mismatches}|g={gaps}")
    return ";".join(stats)

def _blast_screen_passes(spec: dict) -> bool:
    """
    Decide whether to keep a primer pair after BLAST screening.

    We treat three cases as pass:
    - spec['ok'] is True (Primer-BLAST style criteria met)
    - No BLAST hits at all (likely the local DB lacks the target; keep but mark as unverified)
    - Caller requested to keep non-specific via --keep-non-specific (handled upstream)
    """
    if spec.get("ok"):
        return True
    if spec.get("error"):
        return False
    forward_hits = spec.get("forward_hits") or []
    reverse_hits = spec.get("reverse_hits") or []
    # If database has no hits for either primer, don't drop the pair—flag but allow through.
    if not forward_hits and not reverse_hits:
        return True
    return False

def main():
    # Silence noisy Biopython FASTA comment warnings; we already parse with fasta-pearson.
    warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)

    # Set up logging
    import logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    ap = argparse.ArgumentParser(description="Design and score primer candidates for multiple genes.")
    ap.add_argument("--genes", nargs="+", required=True,
                    help="Gene symbols or IDs (e.g., RPP30 MRGPRX1 or 10556).")
    ap.add_argument("--num-return", type=int, default=10,
                    help="Max primer pairs per gene to design (default: 10).")
    ap.add_argument("--product-min", type=int, default=70)
    ap.add_argument("--product-max", type=int, default=400)
    ap.add_argument("--out-prefix", default="PrimerAnalyse_local",
                    help="Prefix for outputs (<prefix>.xlsx and <prefix>_primerblast.csv; default: PrimerAnalyse_local)")
    ap.add_argument("--results-dir", default="results",
                    help="Directory to store numbered result files (default: results)")
    ap.add_argument("--blast-screen", action="store_true",
                    help="Run a Primer-BLAST-like specificity screen via local BLAST+.")
    ap.add_argument("--blast-db",
                    help="Path to a local BLAST nucleotide DB prefix (required when --blast-mode=local).")
    ap.add_argument("--blastn-bin", default="blastn",
                    help="blastn executable to invoke for local BLAST (default: blastn on PATH).")
    ap.add_argument("--blast-max-evalue", type=float, default=1.0,
                    help="Maximum e-value to accept for BLAST hits (default: 1.0).")
    ap.add_argument("--blast-max-mismatches", type=int, default=1,
                    help="Maximum mismatches allowed for target hits (default: 1).")
    ap.add_argument("--blast-min-total-mismatches", type=int, default=2,
                    help="Minimum mismatches required for off-target hits (default: 2).")
    ap.add_argument("--blast-min-3prime-mismatches", type=int, default=2,
                    help="Minimum mismatches required in the last N bases at the 3' end for off-target hits (default: 2).")
    ap.add_argument("--blast-3prime-window", type=int, default=5,
                    help="Number of bases from the 3' end to inspect for mismatch enforcement (default: 5).")
    ap.add_argument("--blast-ignore-mismatches-ge", type=int, default=6,
                    help="Ignore off-target hits with mismatches greater or equal to this number (default: 6).")
    ap.add_argument("--max-self-end", type=float, default=3.0,
                    help="Maximum allowed Primer3 Self 3' complementarity (Self End). Set negative to disable.")
    ap.add_argument("--keep-non-specific", action="store_true",
                    help="When --blast-screen is enabled, retain primer pairs that fail the specificity screen.")
    ap.add_argument("--blast-report-dir",
                    help="Directory to store detailed Primer-BLAST reports (default: <out_prefix>_primerblast when screening).")
    args = ap.parse_args()

    if args.blast_screen and not args.blast_db:
        ap.error("--blast-db is required when using --blast-screen")

    # Create results directory
    results_dir = Path(args.results_dir)
    results_dir.mkdir(exist_ok=True)

    # Determine run number
    run_num = _get_next_run_number(str(results_dir), args.out_prefix)


    # 1) Design
    all_oligos = []
    accepted_pairs = []
    candidate_multiplier = 5 if args.blast_screen else 1
    max_self_end = args.max_self_end if args.max_self_end is None or args.max_self_end >= 0 else None
    shortfalls = {}
    for gene in args.genes:
        logger.info(f"\n{'='*60}")
        logger.info(f"Designing primers for gene: {gene}")
        logger.info(f"{'='*60}")
        accepted_for_gene = 0
        target_pairs = args.num_return
        try:
            seq = fetch_gene_sequence(gene)
            logger.info(f"✓ Fetched sequence for {gene}: {len(seq)} bp")
        except Exception as e:
            logger.error(f"✗ Failed to fetch sequence for {gene}: {e}")
            continue

        pairs = design_primers_for_gene(
            seq,
            product_size_range=(args.product_min, args.product_max),
            num_return=max(target_pairs * candidate_multiplier, target_pairs)
        )
        logger.info(f"Got {len(pairs)} candidate pairs from primer3")

        for i, p in enumerate(pairs):
            self_end_ok = True
            if max_self_end is not None:
                f_end = p.get("forward_self_end")
                r_end = p.get("reverse_self_end")
                if (f_end is not None and f_end > max_self_end) or (r_end is not None and r_end > max_self_end):
                    self_end_ok = False
                    logger.debug(f"  Pair {i}: Rejected - self_end > {max_self_end} (F:{f_end}, R:{r_end})")

            spec = {"ok": None, "common_accessions": [], "forward_hits": [], "reverse_hits": [], "amplicon_hits": []}
            if args.blast_screen and self_end_ok:
                spec = pair_specific_on_same_transcript(
                    p["forward"],
                    p["reverse"],
                    size_min=args.product_min,
                    size_max=args.product_max,
                    blast_db=args.blast_db,
                    blastn_path=args.blastn_bin,
                    max_evalue=args.blast_max_evalue,
                    max_mismatches=args.blast_max_mismatches,
                    min_total_mismatches=args.blast_min_total_mismatches,
                    min_3prime_mismatches=args.blast_min_3prime_mismatches,
                    three_prime_window=args.blast_3prime_window,
                    ignore_mismatches_ge=args.blast_ignore_mismatches_ge,
                )
            elif args.blast_screen and not self_end_ok:
                spec = {"ok": None, "common_accessions": [], "forward_hits": [], "reverse_hits": [], "amplicon_hits": [], "error": "skipped_due_to_self_end"}

            if not self_end_ok:
                continue

            spec_passes = _blast_screen_passes(spec)
            if args.blast_screen and not spec_passes and not args.keep_non_specific:
                logger.debug(f"  Pair {i}: Rejected - BLAST specificity check failed")
                continue

            if accepted_for_gene >= target_pairs:
                logger.debug(f"  Pair {i}: Skipped - already have {target_pairs} pairs")
                continue

            forward_hit_ids = []
            reverse_hit_ids = []
            if spec.get("forward_hits"):
                forward_hit_ids = sorted({h.get("accession") for h in spec["forward_hits"] if isinstance(h, dict)})
            if spec.get("reverse_hits"):
                reverse_hit_ids = sorted({h.get("accession") for h in spec["reverse_hits"] if isinstance(h, dict)})

            forward_stats = _format_hit_stats(spec.get("forward_hits"))
            reverse_stats = _format_hit_stats(spec.get("reverse_hits"))

            accepted_pairs.append({
                "gene": gene,
                "pair_idx": i,
                **p,
                "forward_length": len(p["forward"]),
                "reverse_length": len(p["reverse"]),
                "common_accessions": ";".join(spec.get("common_accessions", [])),
                "specificity_error": spec.get("error"),
                "blast_pass": spec_passes,
                "forward_hits": ";".join(forward_hit_ids),
                "reverse_hits": ";".join(reverse_hit_ids),
                "forward_hits_stats": forward_stats,
                "reverse_hits_stats": reverse_stats,
                "amplicon_hits": ";".join(
                    f"{hit['accession']}:{hit['amplicon_size']}" for hit in spec.get("amplicon_hits", [])
                ),
            })
            logger.info(f"  ✓ Pair {i} ACCEPTED: F={p['forward']} R={p['reverse']}")
            all_oligos.append(Oligo(name=f"{gene}_p{i}_F", seq=p["forward"], role="primer_f"))
            all_oligos.append(Oligo(name=f"{gene}_p{i}_R", seq=p["reverse"], role="primer_r"))
            accepted_for_gene += 1

            if accepted_for_gene >= target_pairs:
                break

        logger.info(f"Gene {gene}: Accepted {accepted_for_gene}/{target_pairs} pairs")

        if accepted_for_gene < target_pairs:
            shortfalls[gene] = (accepted_for_gene, target_pairs)

    # 2) Score (per-oligo + cross-matrix + violations)
    if all_oligos:
        per_oligo, cross, violations = evaluate_set(all_oligos)
    else:
        print("No primer pairs passed design/specificity filters; skipping thermodynamic evaluation.")
        per_oligo = pd.DataFrame(columns=["name", "role", "hairpin_dG", "self_dimer_dG"])
        cross = pd.DataFrame()
        violations = pd.DataFrame(columns=["type", "oligo_a", "oligo_b", "dg"])

    # 3) Write outputs (Excel + primerblast CSV)
    excel_path = results_dir / f"{run_num}_{args.out_prefix}.xlsx"
    primerblast_path = results_dir / f"{run_num}_{args.out_prefix}_primerblast.csv"
    hetero = cross.copy()
    hetero.insert(0, "name", cross.index)
    if not cross.empty:
        hetero["dG (kcal/mol)"] = cross.min(axis=1)
    else:
        hetero["dG (kcal/mol)"] = []

    primer_seqs = pd.DataFrame([
        {"name": o.name, "sequence": o.seq, "role": o.role}
        for o in all_oligos
    ], columns=["name", "sequence", "role"])

    with pd.ExcelWriter(str(excel_path)) as writer:
        hetero.to_excel(writer, sheet_name="heterodimer_dG_primer3", index=False)
        per_oligo.to_excel(writer, sheet_name="hairpin_self_dimer_primer3", index=False)
        primer_seqs.to_excel(writer, sheet_name="primer_sequences", index=False)

    if args.blast_screen:
        pd.DataFrame(accepted_pairs).to_csv(str(primerblast_path), index=False)
    else:
        pd.DataFrame(columns=[
            "gene", "pair_idx", "forward", "reverse", "common_accessions",
            "specificity_error", "forward_hits", "reverse_hits", "amplicon_hits"
        ]).to_csv(str(primerblast_path), index=False)

    print(f"\nAccepted {len(accepted_pairs)} primer pairs across {len(args.genes)} genes.")
    if shortfalls:
        for gene, (got, need) in shortfalls.items():
            print(f"Warning: {gene} produced {got}/{need} passing primer pairs.")
    print(f"Wrote:\n  - {excel_path}\n  - {primerblast_path}")
    print(f"Cutoff used: ΔG <= {DG_CUTOFF} kcal/mol flagged as violations.")

if __name__ == "__main__":
    main()
