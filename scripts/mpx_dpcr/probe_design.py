from mpx_dpcr.evaluate import Oligo, evaluate_set
from Bio.Seq import Seq
import primer3

def design_probe_for_amplicon(template: str, forward: str, reverse: str,
                              num_return: int = 3):
    """
    Find and design up to num_return probe sequences inside the amplicon
    defined by forward and reverse primers. Works even if primers bind to
    reverse complement strand.
    """
    template = template.upper()
    forward = forward.upper()
    reverse = reverse.upper()

    # Try to find primers in the template or its reverse complement
    start = template.find(forward)
    end = template.find(reverse)

    if start == -1 or end == -1:
        rc = str(Seq(template).reverse_complement())
        start = rc.find(forward)
        end = rc.find(reverse)
        if start == -1 or end == -1:
            raise ValueError("Primer sequences not found in template or reverse complement.")
        template = rc  # switch to reverse-complement template

    # If reverse primer comes before forward (opposite orientation)
    if end < start:
        start, end = end, start

    amplicon = template[start + len(forward): end]
    candidates = []
    for i in range(0, len(amplicon) - 20):
        for length in range(20, 31):
            if i + length > len(amplicon):
                break
            seq = amplicon[i:i + length]
            if seq[0] == "G":
                continue
            tm = primer3.calc_tm(seq)
            if 66 <= tm <= 72:
                ol = Oligo(name=f"probe_{i}_{length}", seq=seq, role="probe")
                per, _, v = evaluate_set([ol])
                if v.empty and len(per):
                    candidates.append({
                        "seq": seq,
                        "Tm": tm,
                        "hairpin_dG": per.iloc[0]["hairpin_dG"],
                        "self_dimer_dG": per.iloc[0]["self_dimer_dG"],
                    })

    ranked = sorted(candidates, key=lambda x: abs(x["Tm"] - 69))[:num_return]
    return ranked
