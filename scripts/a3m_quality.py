import os
import math
import glob
import shutil
import subprocess
from collections import Counter
from typing import List, Tuple, Optional, Dict

# =========================
# CONFIG
# =========================
INPUT_DIR = "/content/a3mfiles"         # <-- change to desired directory
OUTPUT_ASM = "/content/tst.tsv"   # output TSV-like

ID_PERCENT = 62
ID_THRESH = ID_PERCENT / 100.0

INCLUDE_DOT_AS_GAP = True

GAP_CHARS = set(["-"] + (["."] if INCLUDE_DOT_AS_GAP else []))

# =========================
# IO + NORMALIZATION
# =========================
def read_fasta_like(path: str) -> Tuple[List[str], List[str]]:
    headers, seqs = [], []
    h, buf = None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    headers.append(h)
                    seqs.append("".join(buf))
                h = line[1:].strip()
                buf = []
            else:
                buf.append(line)
        if h is not None:
            headers.append(h)
            seqs.append("".join(buf))
    return headers, seqs

def strip_a3m_insertions(seq: str) -> str:
    # Remove lowercase letters (insertions in A3M)
    return "".join(ch for ch in seq if not ch.islower())

def normalize_msa_to_query(headers: List[str], seqs: List[str]) -> Tuple[str, List[str]]:
    """
    Returns:
      target_name (first token of query header),
      msa_norm: list of sequences after stripping insertions and pad/trunc to query length
    """
    if not headers or not seqs:
        return "UNKNOWN", []

    seqs2 = [strip_a3m_insertions(s) for s in seqs]
    q = seqs2[0]
    L = len(q)
    if L == 0:
        return headers[0].split()[0] if headers[0] else "UNKNOWN", []

    msa_norm = []
    for s in seqs2:
        if len(s) < L:
            s = s + "-" * (L - len(s))
        elif len(s) > L:
            s = s[:L]
        msa_norm.append(s)

    target = headers[0].split()[0] if headers[0] else "UNKNOWN"
    return target, msa_norm

# =========================
# RAW-MSA METRICS (ONLY)
# =========================
def col_entropy(col: List[str]) -> float:
    residues = [c for c in col if c not in GAP_CHARS]
    if not residues:
        return float("nan")
    counts = Counter(residues)
    n = sum(counts.values())
    ent = 0.0
    for v in counts.values():
        p = v / n
        ent -= p * math.log2(p)
    return ent

def mean_entropy(msa: List[str]) -> float:
    if not msa:
        return float("nan")
    N, L = len(msa), len(msa[0])
    ents = []
    for i in range(L):
        col = [msa[r][i] for r in range(N)]
        e = col_entropy(col)
        if not math.isnan(e):
            ents.append(e)
    return sum(ents) / len(ents) if ents else float("nan")

def gap_fraction(msa: List[str]) -> float:
    if not msa:
        return float("nan")
    N, L = len(msa), len(msa[0])
    gaps = sum(sum(1 for c in s if c in GAP_CHARS) for s in msa)
    return gaps / (N * L)

def pairwise_identity_gapmasked(a: str, b: str) -> float:
    matches, comps = 0, 0
    for x, y in zip(a, b):
        if x in GAP_CHARS or y in GAP_CHARS:
            continue
        comps += 1
        if x == y:
            matches += 1
    return (matches / comps) if comps else 0.0

def avg_identity_to_query(msa: List[str]) -> float:
    if not msa or len(msa) < 2:
        return float("nan")
    q = msa[0]
    vals = []
    for s in msa[1:]:
        comps = sum(1 for x, y in zip(q, s) if x not in GAP_CHARS and y not in GAP_CHARS)
        if comps == 0:
            continue
        vals.append(pairwise_identity_gapmasked(q, s))
    return sum(vals) / len(vals) if vals else float("nan")

# =========================
# NEff_id62 (count only)
# =========================
def hhfilter_available() -> bool:
    return shutil.which("hhfilter") is not None

def count_headers_in_a3m(path: str) -> int:
    n = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n

def run_hhfilter_count_only(a3m_in: str, id_percent: int = 62) -> Optional[int]:
    """
    Runs hhfilter and returns number of sequences in output (count of headers).
    Does NOT return/propagate the filtered MSA to metric computations.
    """
    if not hhfilter_available():
        return None

    out_path = a3m_in + f".hhfilter_id{id_percent}.tmp.a3m"
    try:
        cmd = ["hhfilter", "-i", a3m_in, "-o", out_path, "-id", str(id_percent)]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0 or not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
            return None
        n_hat = count_headers_in_a3m(out_path)
        return n_hat
    finally:
        # Always cleanup
        try:
            if os.path.exists(out_path):
                os.remove(out_path)
        except Exception:
            pass

def greedy_nr_count(msa: List[str], thresh: float = 0.62) -> int:
    kept = []
    for s in msa:
        redundant = False
        for k in kept:
            if pairwise_identity_gapmasked(s, k) >= thresh:
                redundant = True
                break
        if not redundant:
            kept.append(s)
    return len(kept)

def neff_id62_from_raw_msa_only(raw_msa: List[str], a3m_path: str) -> float:
    """
    NEff_id62 = (N_hat at 62%) / L
    - N_hat comes from hhfilter if available, else greedy approximation.
    - L is query length from raw_msa (already normalized).
    """
    if not raw_msa:
        return float("nan")
    L = len(raw_msa[0])
    if L == 0:
        return float("nan")

    n_hat = run_hhfilter_count_only(a3m_path, id_percent=ID_PERCENT)
    if n_hat is None:
        n_hat = greedy_nr_count(raw_msa, thresh=ID_THRESH)
    return n_hat / L

# =========================
# MAIN
# =========================
def compute_metrics_for_file(a3m_path: str) -> Optional[Dict[str, float]]:
    headers, seqs = read_fasta_like(a3m_path)
    if not headers or not seqs:
        return None

    target, raw_msa = normalize_msa_to_query(headers, seqs)
    if not raw_msa:
        return None

    # IMPORTANT: ALL metrics except NEff are computed on RAW MSA ONLY.
    out = {
        "target": target,
        "NEff_id62": neff_id62_from_raw_msa_only(raw_msa, a3m_path),
        "MeanEntropy": mean_entropy(raw_msa),
        "GapFraction": gap_fraction(raw_msa),
        "AvgIdentityToQuery": avg_identity_to_query(raw_msa),
    }
    return out

def main(input_dir: str, output_path: str):
    paths = sorted(glob.glob(os.path.join(input_dir, "*.a3m")))
    if not paths:
        raise FileNotFoundError(f"No .a3m files found in: {input_dir}")

    rows = []
    for p in paths:
        m = compute_metrics_for_file(p)
        if m is not None:
            rows.append(m)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("target\tNEff_id62\tMeanEntropy\tGapFraction\tAvgIdentityToQuery\n")
        for r in rows:
            f.write(
                f"{r['target']}\t"
                f"{r['NEff_id62']:.6f}\t"
                f"{r['MeanEntropy']:.6f}\t"
                f"{r['GapFraction']:.6f}\t"
                f"{r['AvgIdentityToQuery']:.6f}\n"
            )

    print(f"Done. Wrote {len(rows)} rows to: {output_path}")
    if hhfilter_available():
        print(f"hhfilter detected: NEff_id62 uses hhfilter -id {ID_PERCENT} (count-only).")
    else:
        print("hhfilter not found: NEff_id62 uses greedy 62% identity approximation.")

# Run
main(INPUT_DIR, OUTPUT_ASM)
