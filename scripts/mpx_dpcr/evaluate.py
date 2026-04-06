from dataclasses import dataclass
import pandas as pd
import primer3
import itertools

DG_CUTOFF = -25.0

@dataclass
class Oligo:
    name: str
    seq: str
    role: str

def _hp_dG(seq: str) -> float:
    return primer3.calcHairpin(seq).dg / 1000       # convert from cal to kcal

def _self_dimer_dG(seq: str) -> float:
    return primer3.calcHomodimer(seq).dg / 1000

def _hetero_dimer_dG(a: str, b: str) -> float:
    return primer3.calcHeterodimer(a, b).dg / 1000 

def oligo_properties(oligo: Oligo):
    return {
        "hairpin_dG": _hp_dG(oligo.seq),
        "self_dimer_dG": _self_dimer_dG(oligo.seq)
    }

def cross_matrix(oligos):
    names = [o.name for o in oligos]
    mat = pd.DataFrame(index=names, columns=names, dtype=float)
    for i, oi in enumerate(oligos):
        for j, oj in enumerate(oligos):
            if j < i:
                mat.iat[i, j] = mat.iat[j, i]
                continue
            if i == j:
                mat.iat[i, j] = float("nan")
            else:
                mat.iat[i, j] = _hetero_dimer_dG(oi.seq, oj.seq)
    return mat

def evaluate_set(oligos):
    rows = []
    for o in oligos:
        props = oligo_properties(o)
        rows.append({"name": o.name, "role": o.role, **props})
    per_oligo = pd.DataFrame(rows)
    cross = cross_matrix(oligos)
    viol_rows = []
    for _, r in per_oligo.iterrows():
        if r["hairpin_dG"] <= DG_CUTOFF:
            viol_rows.append({"type": "hairpin", "oligo_a": r["name"], "oligo_b": None, "dg": r["hairpin_dG"]})
        if r["self_dimer_dG"] <= DG_CUTOFF:
            viol_rows.append({"type": "self_dimer", "oligo_a": r["name"], "oligo_b": None, "dg": r["self_dimer_dG"]})
    for i, j in itertools.combinations(range(len(oligos)), 2):
        a, b = oligos[i], oligos[j]
        dg = cross.iat[i, j]
        if dg <= DG_CUTOFF:
            viol_rows.append({"type": "hetero_dimer", "oligo_a": a.name, "oligo_b": b.name, "dg": dg})
    violations = pd.DataFrame(viol_rows)
    return per_oligo, cross, violations

def pass_fail_summary(violations):
    return {
        "passes_cutoff": violations.empty,
        "n_violations": 0 if violations.empty else len(violations),
    }
