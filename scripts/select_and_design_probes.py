# scripts/select_and_design_probes.py
import pandas as pd
from mpx_dpcr.fetch import fetch_gene_sequence
from mpx_dpcr.probe_design import design_probe_for_amplicon

def select_and_design():
    # 1. Load primer candidates (produced earlier)
    df = pd.read_csv("multiplex_candidates.csv")
    genes = sorted(df["gene"].unique())

    chosen_rows = []

    print("\nAvailable genes:")
    print(", ".join(genes))
    print("\nFor each gene, you’ll see the available primer pairs with their product sizes.\n")

    for gene in genes:
        sub = df[df["gene"] == gene].copy()
        print(f"==== {gene} ====")
        for pid in sorted(sub["pair_idx"].unique()):
            subset = sub[sub["pair_idx"] == pid]
            f = subset.loc[subset["role"] == "F", "sequence"].values[0]
            r = subset.loc[subset["role"] == "R", "sequence"].values[0]
            print(f"  pair_idx {pid} → F: {f[:8]}... R: {r[:8]}...")
        sel = input(f"Enter pair_idx values for {gene} to keep (comma-separated, e.g., 0,2): ")
        if not sel.strip():
            continue
        keep = [int(x.strip()) for x in sel.split(",") if x.strip().isdigit()]

        # 2. Fetch the full gene sequence
        template = fetch_gene_sequence(gene)

        for pid in keep:
            subset = sub[sub["pair_idx"] == pid]
            f = subset.loc[subset["role"] == "F", "sequence"].values[0]
            r = subset.loc[subset["role"] == "R", "sequence"].values[0]
            print(f"Designing probes for {gene} pair {pid}...")
            probes = design_probe_for_amplicon(template, f, r, num_return=3)
            for i, p in enumerate(probes, 1):
                chosen_rows.append({
                    "gene": gene,
                    "pair_idx": pid,
                    "forward": f,
                    "reverse": r,
                    "probe_idx": i,
                    "probe_seq": p["seq"],
                    "probe_Tm": round(p["Tm"], 2),
                    "probe_hairpin_dG": round(p["hairpin_dG"], 2),
                    "probe_self_dimer_dG": round(p["self_dimer_dG"], 2)
                })

    # 3. Save everything
    if chosen_rows:
        out = pd.DataFrame(chosen_rows)
        out.to_csv("selected_with_probes.csv", index=False)
        print(f"\n✅ Wrote {len(out)} probe candidates to selected_with_probes.csv")
    else:
        print("No probes designed (no selections made).")

if __name__ == "__main__":
    select_and_design()
