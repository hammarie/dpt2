"""Primer-BLAST style specificity checking using local BLAST+."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class BlastHit:
    accession: str
    evalue: float
    pident: float
    length: int
    mismatches: int
    strand: Optional[str]
    start: int
    end: int
    qseq: str
    sseq: str


def _resolve_strand(strand_value, start: int, end: int) -> Optional[str]:
    if strand_value in (None, ""):
        pass
    elif isinstance(strand_value, str):
        val = strand_value.lower()
        if val in {"plus", "+"}:
            return "plus"
        if val in {"minus", "-"}:
            return "minus"
    else:
        try:
            strand_int = int(strand_value)
            if strand_int > 0:
                return "plus"
            if strand_int < 0:
                return "minus"
        except Exception:  # pragma: no cover
            pass
    if start <= end:
        return "plus"
    if start > end:
        return "minus"
    return None


def _mismatches_last_n(qseq: str, sseq: str, n: int) -> int:
    """Count mismatches within the last `n` query bases (3' end)."""
    qseq = qseq.upper()
    sseq = sseq.upper()
    count = 0
    considered = 0
    idx = len(qseq) - 1
    while idx >= 0 and considered < n:
        q_base = qseq[idx]
        s_base = sseq[idx] if idx < len(sseq) else "-"
        if q_base != "-":
            considered += 1
            if s_base == "-" or q_base != s_base:
                count += 1
        idx -= 1
    return count


def blast_primer_local_db(
    primer_seq: str,
    *,
    db_path: str,
    blastn_path: str = "blastn",
    task: str = "blastn-short",
    max_evalue: float = 1.0,
    max_mismatches: int = 1,
) -> List[BlastHit]:
    """Run blastn locally against a nucleotide DB and return detailed hits."""
    fasta = f">primer\n{primer_seq}\n"
    cmd = [
        blastn_path,
        "-task",
        task,
        "-db",
        db_path,
        "-outfmt",
        "6 sacc evalue pident length mismatch gaps sstrand sstart send qseq sseq",
        "-query",
        "-",
    ]
    proc = subprocess.run(
        cmd,
        input=fasta,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"blastn failed: {proc.stderr.strip() or proc.stdout.strip()}")

    hits: List[BlastHit] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 9:
            continue
        acc, evalue, pident, length, mismatch, gaps, strand, sstart, send, qseq, sseq = parts
        length = int(float(length))
        mismatches = int(float(mismatch))
        evalue = float(evalue)
        if evalue > max_evalue or mismatches > max_mismatches:
            continue
        start = int(float(sstart))
        end = int(float(send))
        strand_label = _resolve_strand(strand, start, end)
        hits.append(
            BlastHit(
                accession=acc,
                evalue=evalue,
                pident=float(pident),
                length=length,
                mismatches=mismatches,
                strand=strand_label,
                start=start,
                end=end,
                qseq=qseq,
                sseq=sseq,
            )
        )
    return hits


def _passes_specificity_stringency(
    hits: List[BlastHit],
    allowed_accessions: set,
    min_total_mismatches: int,
    min_3prime_mismatches: int,
    three_prime_window: int,
    ignore_mismatches_ge: int,
) -> bool:
    for hit in hits:
        if hit.accession in allowed_accessions:
            continue
        total = hit.mismatches
        if ignore_mismatches_ge is not None and total >= ignore_mismatches_ge:
            continue
        end_mism = _mismatches_last_n(hit.qseq, hit.sseq, three_prime_window)
        if total < min_total_mismatches or end_mism < min_3prime_mismatches:
            return False
    return True


def pair_specific_on_same_transcript(
    forward: str,
    reverse: str,
    *,
    size_min: int = 60,
    size_max: int = 120,
    blast_db: Optional[str] = None,
    blastn_path: str = "blastn",
    max_evalue: float = 1.0,
    max_mismatches: int = 1,
    min_total_mismatches: int = 2,
    min_3prime_mismatches: int = 2,
    three_prime_window: int = 5,
    ignore_mismatches_ge: int = 6,
) -> Dict[str, object]:
    """
    Approximate Primer-BLAST's pair-level specificity check.

    A pair passes iff both primers hit the same accession, on opposite strands,
    and the implied amplicon size lies within [size_min, size_max].
    """
    if not blast_db:
        raise ValueError("blast_db path required for BLAST specificity screening")
    try:
        forward_hits = blast_primer_local_db(
            forward,
            db_path=blast_db,
            blastn_path=blastn_path,
            max_evalue=max_evalue,
            max_mismatches=max_mismatches,
        )
        reverse_hits = blast_primer_local_db(
            reverse,
            db_path=blast_db,
            blastn_path=blastn_path,
            max_evalue=max_evalue,
            max_mismatches=max_mismatches,
        )
    except Exception as exc:  # pragma: no cover - defensive
        return {"ok": False, "common_accessions": [], "error": str(exc)}

    amp_candidates = []
    common_accessions = set()
    for fh in forward_hits:
        for rh in reverse_hits:
            if fh.accession != rh.accession:
                continue
            if fh.strand != "plus" or rh.strand != "minus":
                continue
            f_start = min(fh.start, fh.end)
            f_end = max(fh.start, fh.end)
            r_start = min(rh.start, rh.end)
            r_end = max(rh.start, rh.end)
            if r_end <= f_start:
                continue
            amplicon = r_end - f_start + 1
            if not (size_min <= amplicon <= size_max):
                continue
            common_accessions.add(fh.accession)
            amp_candidates.append(
                {
                    "accession": fh.accession,
                    "forward_start": f_start,
                    "forward_end": f_end,
                    "reverse_start": r_start,
                    "reverse_end": r_end,
                    "amplicon_size": amplicon,
                    "forward_evalue": fh.evalue,
                    "reverse_evalue": rh.evalue,
                }
            )

    allowed = set(common_accessions)
    specificity_ok = False
    if amp_candidates and allowed:
        specificity_ok = _passes_specificity_stringency(
            forward_hits,
            allowed,
            min_total_mismatches,
            min_3prime_mismatches,
            three_prime_window,
            ignore_mismatches_ge,
        ) and _passes_specificity_stringency(
            reverse_hits,
            allowed,
            min_total_mismatches,
            min_3prime_mismatches,
            three_prime_window,
            ignore_mismatches_ge,
        )

    ok = bool(amp_candidates) and specificity_ok

    return {
        "ok": ok,
        "common_accessions": sorted(common_accessions),
        "forward_hits": [fh.__dict__ for fh in forward_hits],
        "reverse_hits": [rh.__dict__ for rh in reverse_hits],
        "amplicon_hits": amp_candidates,
        "size_min": size_min,
        "size_max": size_max,
        "max_evalue": max_evalue,
        "max_mismatches": max_mismatches,
        "blast_db": blast_db,
        "min_total_mismatches": min_total_mismatches,
        "min_3prime_mismatches": min_3prime_mismatches,
        "three_prime_window": three_prime_window,
        "ignore_mismatches_ge": ignore_mismatches_ge,
        "specificity_stringency_ok": specificity_ok,
    }
