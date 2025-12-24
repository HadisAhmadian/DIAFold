#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clean monomeric protein mmCIFs for structure-prediction benchmarks.

Assumptions (enforced upstream by your RCSB query):
  - Entry contains a single protein chain (monomeric, protein-only).
  - No DNA/RNA polymers.
  - Any major "investigated" ligands or covalent linkers already filtered out.

What this script does:
  1. Loads each .cif file and its mmCIF dictionary.
  2. Selects model 0 and the first chain in that model as the protein chain.
  3. Uses _pdbx_poly_seq_scheme to:
       - Get the polymer sequence positions for that chain.
       - Determine which positions have coordinates (coverage).
       - Compute max internal gap on the polymer sequence.
  4. Rejects entries with:
       - coverage < 0.8
       - max internal gap > 15
  5. Cleans the structure:
       - Keeps only the chosen chain.
       - Within that chain, keeps only residues that belong to the polymer
         sequence (as defined by _pdbx_poly_seq_scheme).
       - This removes waters, ions, non-polymer ligands, etc.
       - PTMs that are part of the polymer sequence are kept automatically.
  6. Writes cleaned .cif files and a simple TSV log.

The log will include 'method:main' or 'method:fallback' for
coverage/gap-based decisions so you can see which path was used.
"""

import os
import glob
import argparse

from Bio.PDB import MMCIFParser, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # IMPORTANT: this is the correct import

COVERAGE_CUTOFF = 0.8
MAX_INTERNAL_GAP = 15


# ---------------------------------------------------------------------------

def _to_list(x):
    """Ensure mmCIF field is a list (MMCIF2Dict returns str for singletons)."""
    if isinstance(x, list):
        return x
    return [x]


def get_poly_seq_scheme(cif_dict, chain_id):
    """
    Extract polymer sequence information for a chain from _pdbx_poly_seq_scheme.

    Returns:
      seq_ids: sorted list of unique polymer positions (integers)
      auth_map: dict mapping (auth_seq_num, ins_code) -> set(seq_id)

    Only rows where:
      - pdb_strand_id == chain_id
      - seq_id is an integer

    are used.
    """
    try:
        seq_ids_raw = _to_list(cif_dict["_pdbx_poly_seq_scheme.seq_id"])
        strand_ids = _to_list(cif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"])
        auth_nums_raw = _to_list(
            cif_dict.get("_pdbx_poly_seq_scheme.auth_seq_num",
                         ["?"] * len(seq_ids_raw))
        )
        ins_codes_raw = _to_list(
            cif_dict.get("_pdbx_poly_seq_scheme.pdb_ins_code",
                         ["?"] * len(seq_ids_raw))
        )
    except KeyError:
        return None, None

    auth_map = {}
    seq_id_set = set()

    for s_id, st_id, a_num, ins_code in zip(
        seq_ids_raw, strand_ids, auth_nums_raw, ins_codes_raw
    ):
        if st_id.strip() != chain_id:
            continue

        # polymer position
        try:
            seq_id = int(s_id)
        except Exception:
            continue

        seq_id_set.add(seq_id)

        # author residue number (may be ".", "?", etc.)
        try:
            auth_num = int(a_num)
        except Exception:
            continue

        icode = "" if ins_code in (".", "?", None) else ins_code.strip()
        key = (auth_num, icode)
        auth_map.setdefault(key, set()).add(seq_id)

    if not seq_id_set:
        return None, None

    seq_ids = sorted(seq_id_set)
    return seq_ids, auth_map


def compute_coverage_and_gap(chain, seq_ids, auth_map):
    """
    Compute coverage and max internal gap on the polymer sequence.

    Args:
      chain: Biopython Chain object (chosen protein chain)
      seq_ids: sorted list of polymer positions for this chain
      auth_map: mapping (auth_seq_num, icode) -> set(seq_id)

    Returns:
      coverage (float), max_internal_gap (int)
    """
    if not seq_ids:
        return 0.0, 0

    observed_seq_ids = set()

    for res in chain.get_residues():
        het_flag, resseq, icode = res.id
        icode_str = "" if icode == " " else str(icode)
        key = (resseq, icode_str)
        if key in auth_map:
            observed_seq_ids.update(auth_map[key])

    if not observed_seq_ids:
        return 0.0, 0

    # coverage
    coverage = len(observed_seq_ids) / len(seq_ids)

    # max internal gap
    obs_sorted = sorted(observed_seq_ids)
    if len(obs_sorted) < 2:
        return coverage, 0

    max_gap = 0
    prev = obs_sorted[0]
    for pos in obs_sorted[1:]:
        gap = pos - prev - 1
        if gap > max_gap:
            max_gap = gap
        prev = pos

    return coverage, max_gap


def fallback_coverage_gap(chain):
    """
    Fallback estimate if _pdbx_poly_seq_scheme is missing or unusable.

    Uses simple author residue numbers (resseq) to approximate coverage
    and internal gaps. This is less accurate but should rarely be used
    for modern mmCIFs.
    """
    resseqs = [res.id[1] for res in chain.get_residues()]
    if not resseqs:
        return 0.0, 0

    resseqs_sorted = sorted(set(resseqs))
    min_r, max_r = resseqs_sorted[0], resseqs_sorted[-1]
    length_est = max_r - min_r + 1
    coverage = len(resseqs_sorted) / length_est if length_est > 0 else 0.0

    max_gap = 0
    prev = resseqs_sorted[0]
    for r in resseqs_sorted[1:]:
        gap = r - prev - 1
        if gap > max_gap:
            max_gap = gap
        prev = r

    return coverage, max_gap


def clean_chain_using_polymer_mapping(chain, auth_map):
    """
    Remove all residues from the chain that are NOT part of the polymer
    sequence defined by auth_map.

    Any residue whose (auth_seq_num, icode) appears in auth_map is kept.
    Everything else (waters, ions, ligands, extra junk) is removed.
    """
    for res in list(chain.get_residues()):
        het_flag, resseq, icode = res.id
        icode_str = "" if icode == " " else str(icode)
        key = (resseq, icode_str)
        if key not in auth_map:
            chain.detach_child(res.id)


# ---------------------------------------------------------------------------

def process_cif_file(cif_path, out_dir, log_fh, parser, io):
    basename = os.path.basename(cif_path)

    # Parse coordinates
    try:
        structure = parser.get_structure("prot", cif_path)
    except Exception as e:
        msg = f"{basename}\tREJECT\tparse_error:{e}\n"
        log_fh.write(msg)
        print(msg.strip())
        return

    # Parse mmCIF dictionary
    try:
        cif_dict = MMCIF2Dict(cif_path)
    except Exception as e:
        cif_dict = None
        print(f"{basename}: warning: MMCIF2Dict failed ({e}). "
              f"Using fallback coverage/gap.")

    # Get model 0
    try:
        model = structure[0]
    except Exception as e:
        msg = f"{basename}\tREJECT\tno_model_0:{e}\n"
        log_fh.write(msg)
        print(msg.strip())
        return

    # Assume single protein chain (enforced by RCSB filter); take first chain
    chains = list(model.get_chains())
    if not chains:
        msg = f"{basename}\tREJECT\tno_chains_found\n"
        log_fh.write(msg)
        print(msg.strip())
        return

    chain = chains[0]
    chain_id = chain.id

    # 1) Coverage & internal gap, track method used
    if cif_dict is not None:
        seq_ids, auth_map = get_poly_seq_scheme(cif_dict, chain_id)
    else:
        seq_ids, auth_map = None, None

    if seq_ids is None or auth_map is None:
        method = "fallback"
        print(f"{basename}: WARNING – using fallback coverage/gap "
              f"(author numbering only)")
        coverage, max_gap = fallback_coverage_gap(chain)
    else:
        method = "main"
        coverage, max_gap = compute_coverage_and_gap(chain, seq_ids, auth_map)

    # Apply thresholds
    if coverage < COVERAGE_CUTOFF:
        msg = (f"{basename}\tREJECT\tlow_coverage:{coverage:.3f}"
               f"\tmethod:{method}\n")
        log_fh.write(msg)
        print(msg.strip())
        return

    if max_gap > MAX_INTERNAL_GAP:
        msg = (f"{basename}\tREJECT\tbig_internal_gap:{max_gap}"
               f"\tmethod:{method}\n")
        log_fh.write(msg)
        print(msg.strip())
        return

    # 2) Clean structure:
    #    - keep only chosen chain
    #    - keep only polymer residues (according to auth_map if available)
    for other_chain in list(model.get_chains()):
        if other_chain.id != chain_id:
            model.detach_child(other_chain.id)

    if auth_map is not None:
        clean_chain_using_polymer_mapping(chain, auth_map)
    else:
        # Fallback: remove all HET residues (waters, ions, ligands)
        for res in list(chain.get_residues()):
            het_flag, resseq, icode = res.id
            if het_flag.strip() != "":
                chain.detach_child(res.id)

    # Save cleaned structure
    out_path = os.path.join(out_dir, basename)
    io.set_structure(structure)
    io.save(out_path)

    msg = (
        f"{basename}\tACCEPT\tok"
        f"\tchain:{chain_id}"
        f"\tcoverage:{coverage:.3f}"
        f"\tmax_gap:{max_gap}"
        f"\tmethod:{method}\n"
    )
    log_fh.write(msg)
    print(msg.strip())


# ---------------------------------------------------------------------------

def main():
    cli = argparse.ArgumentParser(
        description="Clean monomeric protein mmCIFs for structure-prediction benchmarks."
    )
    cli.add_argument(
        "-i", "--input_dir", required=True,
        help="Input directory containing .cif files"
    )
    cli.add_argument(
        "-o", "--output_dir", required=True,
        help="Output directory for cleaned .cif files and log"
    )
    args = cli.parse_args()

    in_dir = args.input_dir
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, "clean_log.txt")

    parser = MMCIFParser(QUIET=True)
    io = MMCIFIO()

    cif_files = sorted(glob.glob(os.path.join(in_dir, "*.cif")))
    if not cif_files:
        print("No .cif files found in input directory:", in_dir)
        return

    with open(log_path, "w", encoding="utf-8") as log_fh:
        log_fh.write("# file_name\tSTATUS\treason_or_metrics\n")
        for cif_path in cif_files:
            process_cif_file(cif_path, out_dir, log_fh, parser, io)

    print("Done. Log written to:", log_path)


if __name__ == "__main__":
    main()
