"""Calculate LMU benchmark ΔG metrics for predefined primer panels."""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path
from typing import Dict, List

import pandas as pd

from mpx_dpcr.evaluate import Oligo, evaluate_set

try:  # Optional dependency for alternative thermodynamics
    from nupack import Complex, Model, Strand, complex_analysis

    NUPACK_AVAILABLE = True
except Exception:  # pragma: no cover - optional
    NUPACK_AVAILABLE = False


PRIMER_SEQUENCES: Dict[str, str] = {
    "IL12A_PP5_FW": "AAGACCTCTTTTATGATGGCCC",
    "IL12A_PP5_RV": "TGGCACAGTCTCACTGTTGA",
    "IL12A_PP10_FW": "GCCGTCAGCAACATGCTC",
    "IL12A_PP10_RV": "CTACTAAGGCACAGGGCCATCATAA",
    "Il10_PP1_FW": "CAAGCTGAGAACCAAGACCCA",
    "IL10_PP1_FW": "CAAGCTGAGAACCAAGACCCA",
    "IL10_PP_RV": "CACAGGGAAGAAATCGATGACAGC",
    "IL10_PP2_FW": "GCTGAGAACCAAGACCCAGA",
    "IL_10_PP2_RV": "ATGCCTTTCTCTTGGAGCTTAT",
    "IL1B_PP1_FW": "AGCTTGGTGATGTCTGGTCC",
    "IL1B_PP1_RV": "TGGAGAACACCACTTGTTGC",
    "IL1B_PP1_Rv": "TGGAGAACACCACTTGTTGC",
    "IL1B_PP10_FW": "ATCTGTACCTGTCCTGCGTG",
    "IL1B_PP10_RV": "TTTTTGGGATCTACACTCTCCAGC",
    "IL-6 FW": "TGGCAGAAAACAACCTGAACC",
    "IL-6RV": "TTTCACCAGGCAAGTCTCCTCAT",
    "IL8_FW": "ACCACCGGAAGGAACCATCT",
    "IL8_RV": "TGGCAAAACTGCACCTTCACA",
    "TNFa_FW": "AAAACAACCCTCAGACGCCA",
    "TNFa_RV": "TCCTTTCCAGGGGAGAGAGG",
    "B2M_FW": "AGATGAGTATGCCTGCCGTG",
    "B2M_RV": "ACCTCCATGATGCTGCTTACA",
    "GDPH_FW": "CCACATCGCTCAGACACCAT",  # GAPDH alias
    "GDPH_RV": "TGAAGGGGTCATTGATGGCAA",
}

COLUMN_NAMES: List[str] = [
    "Il10_PP1_FW",
    "IL10_PP_RV",
    "IL10_PP2_FW",
    "IL_10_PP2_RV",
    "IL1B_PP1_FW",
    "IL1B_PP1_RV",
    "IL1B_PP10_FW",
    "IL1B_PP10_RV",
    "IL-6 FW",
    "IL-6RV",
    "IL8_FW",
    "IL8_RV",
    "TNFa_FW",
    "TNFa_RV",
    "B2M_FW",
    "B2M_RV",
    "GDPH_FW",
    "GDPH_RV",
]

ROW_NAMES: List[str] = [
    "IL12A_PP5_FW",
    "IL12A_PP5_RV",
    "IL12A_PP10_FW",
    "IL12A_PP10_RV",
    "Il10_PP1_FW",
    "IL10_PP_RV",
    "IL10_PP2_FW",
    "IL_10_PP2_RV",
    "IL1B_PP1_FW",
    "IL1B_PP1_RV",
    "IL1B_PP10_FW",
    "IL1B_PP10_RV",
    "IL-6 FW",
    "IL-6RV",
    "IL8_FW",
    "IL8_RV",
    "TNFa_FW",
    "TNFa_RV",
]


UNION_ORDER: List[str] = []
for _name in ROW_NAMES + COLUMN_NAMES:
    if _name not in UNION_ORDER:
        UNION_ORDER.append(_name)


def resolve_sequence(label: str) -> str:
    try:
        seq = PRIMER_SEQUENCES[label]
    except KeyError:
        raise KeyError(f"No sequence available for primer '{label}'")
    return seq.replace(" ", "").upper()


def build_oligos() -> List[Oligo]:
    oligos: List[Oligo] = []
    missing = []
    for label in UNION_ORDER:
        try:
            seq = resolve_sequence(label)
        except KeyError:
            missing.append(label)
            continue
        oligos.append(Oligo(name=label, seq=seq, role="primer"))
    if missing:
        raise ValueError(f"Sequences missing for: {', '.join(missing)}")
    return oligos


def highlight_threshold(val, thresh: float):
    if pd.isna(val):
        return ""
    if abs(val) >= thresh:
        return "background-color:#f8d7da"
    return ""


def nupack_model():
    if not NUPACK_AVAILABLE:
        return None
    return Model(material="dna", celsius=37)


if NUPACK_AVAILABLE:  # pragma: no cover - optional
    NU_MODEL = nupack_model()

    @lru_cache(maxsize=None)
    def _nu_strand(name: str):
        seq = resolve_sequence(name)
        return Strand(seq, name=name)

    def nupack_energy(strand_names: List[str]) -> float:
        strands = [_nu_strand(n) for n in strand_names]
        complex_obj = Complex(strands=strands, name="+".join(strand_names))
        result = complex_analysis(complexes=[complex_obj], model=NU_MODEL)
        return result[complex_obj].free_energy

    def nupack_hairpin(name: str) -> float:
        return nupack_energy([name])

    def nupack_self_dimer(name: str) -> float:
        return nupack_energy([name, name])

    def nupack_heterodimer(name_a: str, name_b: str) -> float:
        return nupack_energy([name_a, name_b])


def sequences_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        [{"name": n, "sequence": resolve_sequence(n)} for n in UNION_ORDER]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute LMU benchmark ΔG matrices.")
    parser.add_argument("--out-xlsx", default="benchmark_lmu.xlsx",
                        help="Output Excel workbook (default: benchmark_lmu.xlsx)")
    parser.add_argument("--highlight-threshold", type=float, default=9.0,
                        help="Absolute ΔG threshold for highlighting cells red (default: 9.0)")
    args = parser.parse_args()

    oligos = build_oligos()
    per_oligo, cross, _ = evaluate_set(oligos)

    cross = cross.reindex(index=UNION_ORDER, columns=UNION_ORDER)
    subset = cross.loc[ROW_NAMES, COLUMN_NAMES]
    per_subset = per_oligo[per_oligo["name"].isin(UNION_ORDER)].copy()
    seq_df = sequences_dataframe()

    nupack_subset = None
    nupack_per = None
    if NUPACK_AVAILABLE:
        try:
            nupack_subset = pd.DataFrame(
                [[nupack_heterodimer(r, c) for c in COLUMN_NAMES] for r in ROW_NAMES],
                index=ROW_NAMES,
                columns=COLUMN_NAMES,
            )
            nupack_per = pd.DataFrame(
                [
                    {
                        "name": name,
                        "hairpin_dG": nupack_hairpin(name),
                        "self_dimer_dG": nupack_self_dimer(name),
                    }
                    for name in UNION_ORDER
                ]
            )
        except Exception as exc:  # pragma: no cover - optional
            print(f"Warning: NUPACK calculations failed: {exc}")
            nupack_subset = None
            nupack_per = None
    else:
        print("NUPACK python package not installed; NUPACK sheet will contain a notice.")

    out_path = Path(args.out_xlsx)
    with pd.ExcelWriter(out_path) as writer:
        subset.style.format("{:.3f}").applymap(
            lambda v: highlight_threshold(v, args.highlight_threshold)
        ).to_excel(writer, sheet_name="heterodimer_dG_primer3")
        per_subset.to_excel(writer, sheet_name="hairpin_self_dimer_primer3", index=False)
        seq_df.to_excel(writer, sheet_name="primer_sequences", index=False)

        if nupack_subset is not None and nupack_per is not None:
            nupack_subset.style.format("{:.3f}").applymap(
                lambda v: highlight_threshold(v, args.highlight_threshold)
            ).to_excel(writer, sheet_name="heterodimer_dG_nupack")
            nupack_per.to_excel(writer, sheet_name="hairpin_self_dimer_nupack", index=False)
        else:
            pd.DataFrame({
                "message": [
                    "NUPACK calculations unavailable. Install the 'nupack' package to enable this sheet."
                ]
            }).to_excel(writer, sheet_name="heterodimer_dG_nupack", index=False)

    print(f"Wrote LMU benchmark workbook to {out_path}")


if __name__ == "__main__":
    main()
