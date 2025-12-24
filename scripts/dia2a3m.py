from Bio import SeqIO
import argparse
from pathlib import Path


def read_fasta_file(file_path):
    print("Reading FASTA:", file_path)
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq).strip()
    print("Done reading FASTA")
    return sequences


def hits_per_query(fname):
    """
    Reads the result file and groups lines by query ID.
    Each value in the dictionary is a list of:
      [seq, cigar, q_start, s_start, header, ...]
    """
    with open(fname, "r") as file1:
        lines = file1.readlines()

    lines = [l.split() for l in lines]
    d = {}

    for l in lines:
        if l[0] not in d:
            d[l[0]] = []
        d[l[0]].append(l[1:])

    return d


def cigar_to_alignment_fold(cigar_string, sequence, query):
    """
    Converts CIGAR string to aligned query and subject sequences.
    Rules (allowDeletion = 0 approach):

    M: match/mismatch -> new column (both take 1 char)
    I: insertion in query -> add char to query / '-' to subject
    D: insertion in subject relative to query -> store lowercase chars in insertions[col]

    Output:
      alignment_q : aligned query string
      alignment_s : aligned subject string
      insertions  : list of lowercase strings per column
    """
    alignment_s = []
    alignment_q = []
    insertions = []

    index = 0
    s_curr = 0
    q_curr = 0
    col = -1

    for char in cigar_string:
        if char.isdigit():
            index = index * 10 + int(char)
        else:
            if char == "M":
                for _ in range(index):
                    if q_curr >= len(query) or s_curr >= len(sequence):
                        break
                    alignment_q.append(query[q_curr])
                    alignment_s.append(sequence[s_curr])
                    insertions.append("")
                    q_curr += 1
                    s_curr += 1
                    col += 1

            elif char == "I":
                for _ in range(index):
                    if q_curr >= len(query):
                        break
                    alignment_q.append(query[q_curr])
                    alignment_s.append("-")
                    insertions.append("")
                    q_curr += 1
                    col += 1

            elif char == "D":
                for _ in range(index):
                    if s_curr >= len(sequence):
                        break
                    if col >= 0:
                        insertions[col] += sequence[s_curr].lower()
                    s_curr += 1

            else:
                raise ValueError("Invalid CIGAR character: " + char)

            index = 0  # reset

    while q_curr < len(query):
        alignment_q.append(query[q_curr])
        alignment_s.append("-")
        insertions.append("")
        q_curr += 1
        col += 1

    return "".join(alignment_q), "".join(alignment_s), insertions


def hits_to_a3m(q_seq, all_hits, outfname, out_dir: Path):
    """
    Creates A3M alignment file for one query.
    """
    q_seq = q_seq.strip()
    target_alignments = []

    for l in all_hits:
        hit_seq = l[0]
        cigar = l[1]
        q_start = int(l[2]) - 1
        s_start = int(l[3]) - 1
        header = l[4]

        q_prefix = q_seq[:q_start]
        q_aln = q_seq[q_start:]

        s_prefix = "-" * len(q_prefix)
        s_sub = hit_seq[s_start:]

        _, s_aln, insertions = cigar_to_alignment_fold(cigar, s_sub, q_aln)

        s_full = s_prefix + s_aln
        ins_full = [""] * len(q_prefix) + insertions

        if len(s_full) < len(q_seq):
            s_full += "-" * (len(q_seq) - len(s_full))
            ins_full += [""] * (len(q_seq) - len(ins_full))
        elif len(s_full) > len(q_seq):
            s_full = s_full[:len(q_seq)]
            ins_full = ins_full[:len(q_seq)]

        target_alignments.append((s_full, ins_full, header))

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{outfname}.a3m"
    print("Writing A3M:", out_path)

    with out_path.open("w", encoding="utf-8") as f:
        f.write(">" + outfname + "\n")
        f.write(q_seq + "\n")
        for s_full, ins_full, header in target_alignments:
            f.write(">" + header + "\n")
            pieces = []
            for ch, ins in zip(s_full, ins_full):
                pieces.append(ch)
                if ins:
                    pieces.append(ins)
            f.write("".join(pieces) + "\n")


# ================== MAIN ================== #

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert FASTA and search hits (CIGAR format) to A3M"
    )
    parser.add_argument("--fasta", required=True, help="FASTA file")
    parser.add_argument("--hits", required=True, help="Hits file with CIGAR")
    parser.add_argument("--out", required=True, help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out).expanduser().resolve()
    if out_dir.exists() and not out_dir.is_dir():
        raise ValueError(f"--out must be a directory path, got a file: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    queries_seq = read_fasta_file(args.fasta)
    dict_hits_per_query = hits_per_query(args.hits)

    for query, all_hits in dict_hits_per_query.items():
        if query not in queries_seq:
            print(f"Warning: Query {query} not found in FASTA, skipping.")
            continue
        hits_to_a3m(queries_seq[query], all_hits, query, out_dir)
