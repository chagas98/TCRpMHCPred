import subprocess, sys, os
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import three_to_one
import argparse
from Bio.PDB.PDBIO import Select

"""
Batch-renumber PDBs by selecting and renumbering the best-aligning chain for MHC-I (single reference).
"""
import os
import argparse
import subprocess
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import three_to_one

# ==== EXISTING FUNCTIONS: DO NOT MODIFY BELOW THIS LINE ====
def _is_backbone_complete(residue, ignore_alts=False):
    """Return True if residue has at least N, CA, C, O (no hydrogens)."""
    needed = {"N", "CA", "C", "O"}
    have = set()
    
    for atom in residue.get_atoms():
        name = atom.get_name()
        if name.startswith("H"):
            continue
        if ignore_alts and atom.get_altloc().strip():
            # skip non-primary altlocs if requested
            continue
        if name in needed and atom.get_vector() is not None:
            have.add(name)

    return needed.issubset(have)

def run_tmalign_with_score(ref_pdb, target_pdb):
    """
    Returns (score, aln1, aln2) from a single TM-align run.
    """
    result = subprocess.run(
        ['TMalign', target_pdb, ref_pdb, '-a'],
        stdout=subprocess.PIPE, text=True
    )
    lines = result.stdout.splitlines()
    score = 0.0
    for L in lines:
        if "TM-score=" in L and "Chain_2" in L:
            try:
                score = float(L.split('= ')[1].split()[0])
            except ValueError:
                score = 0.0
            break
    aln1 = aln2 = ""
    for i, line in enumerate(lines):
        if line.startswith("Aligned"):
            aln1 = lines[i + 7].strip()
            aln2 = lines[i + 9].strip()
            break
    return score, aln1, aln2


def extract_residue_mapping_with_pdb_numbers(ref_pdb, target_pdb, aln1, aln2):
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", ref_pdb)
    tgt_struct = parser.get_structure("tgt", target_pdb)
    # select best chain by sequence overlap
    def get_seq(chain):
        seq = ""
        for r in chain:
            if r.id[0] != ' ': continue
            try:
                seq += three_to_one(r.resname)
            except:
                seq += 'X'
        return seq
    def pick_chain(structure, aln):
        gapless = aln.replace('-', '')
        best, best_score = None, -1
        for chain in structure[0]:
            s = get_seq(chain)
            sc = sum(a == b for a, b in zip(gapless, s))
            if sc > best_score:
                best_score, best = sc, chain.id
        return best
    rc = pick_chain(ref_struct, aln2)
    tc = pick_chain(tgt_struct, aln1)
    print(f"Using chains: ref={rc}, tgt={tc}")
    ref_residues = [r for r in ref_struct[0][rc] if r.id[0] == ' ']
    tgt_residues = [r for r in tgt_struct[0][tc] if r.id[0] == ' ' and _is_backbone_complete(r)]

    print(f"Reference chain {rc} has {len(ref_residues)} standard residues.")
    print(f"Target chain {tc} has {len(tgt_residues)} standard residues.")
    
    # build aligned index lists
    ref_iter = iter(ref_residues)
    tgt_iter = iter(tgt_residues)
    ref_positions = []
    tgt_positions = []
    for a1, a2 in zip(aln1, aln2):
        if a2 != '-':
            r = next(ref_iter)
            print(r)
            ref_positions.append(r.id[1])
        else:
            ref_positions.append(None)
        if a1 != '-':
            t = next(tgt_iter)
            tgt_positions.append(t)
        else:
            tgt_positions.append(None)

    print(f"Alignment length: {len(ref_positions)} residues.")
    non_none = [p for p in ref_positions if p is not None]
    if not non_none:
        return {}
    first_ref = non_none[0]
    last_ref = non_none[-1]
    mapping = {}
    insertion_counters = {}
    total = len(ref_positions)

    # assign mappings
    for idx, (ref_pos, tgt_r) in enumerate(zip(ref_positions, tgt_positions)):
        if tgt_r is None:
            continue
        if ref_pos is not None:
            mapping[(tc, tgt_r.id)] = (ref_pos, rc, ' ')
        else:
            if idx < ref_positions.index(first_ref):
                newnum = first_ref - (ref_positions.index(first_ref) - idx)
                mapping[(tc, tgt_r.id)] = (newnum, rc, ' ')
            elif idx > (total - 1 - ref_positions[::-1].index(last_ref)):
                newnum = last_ref + (idx - (total - 1 - ref_positions[::-1].index(last_ref)))
                mapping[(tc, tgt_r.id)] = (newnum, rc, ' ')
            else:
                prev_idxs = [i for i, p in enumerate(ref_positions[:idx]) if p is not None]
                base = ref_positions[prev_idxs[-1]]
                cnt = insertion_counters.get(base, 0)
                letter = chr(ord('A') + cnt)
                insertion_counters[base] = cnt + 1
                mapping[(tc, tgt_r.id)] = (base, rc, letter)
    print(f"Extracted {len(mapping)} residue mappings.")
    return mapping


def renumber_target_pdb(target_pdb, output_pdb, mapping):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tgt", target_pdb)
    chains_to_renum = set(cid for cid, _ in mapping.keys())
    print(f"Renumbering chains: {chains_to_renum}")

    for model in structure:
        for chain in list(model):
            if chain.id not in chains_to_renum:
                continue
            old_id = chain.id
            residues = [r for r in chain if r.id[0] == ' ']
            max_orig = max(r.id[1] for r in residues)
            max_new = max(num for num, cid, _ in mapping.values() if cid == old_id)
            offset = max_orig + max_new + 1
            new_resids = []
            for r in residues:
                new_resids.append(r.id[1] + offset)
                r.id = (' ', r.id[1] + offset, r.id[2])
            
            keep_set = set()
            for (cid, orig_id), (new_num, new_chain, icode) in mapping.items():
                if cid != old_id:
                    continue
                key = (' ', orig_id[1] + offset, orig_id[2])

                res = next((x for x in residues if x.id == key), None)

                if res is None:
                    # not found -> skip
                    continue
                else:
                    # found -> renumber
                    res.id = (' ', new_num, icode)
                    keep_set.add(res.id[1])

            # drop unmapped standard residues from the chain
            chain.child_list = [
                r for r in chain.child_list
                if (r.id[0] != ' ') or (r.id[1] in keep_set)
            ]

            # update chain ID if all residues map to a single new chain ID
            new_chain_ids = {newc for (_, _), (_, newc, _) in mapping.items() if _ == old_id}
            if len(new_chain_ids) == 1:
                chain.id = new_chain_ids.pop()

    print(f"Renumbered PDB saved to: {output_pdb}")
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)


class ChainSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id


def extract_chain(input_pdb, output_pdb, chain_id):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('TMP', input_pdb)
    io = PDBIO()
    io.set_structure(struct)
    io.save(output_pdb, select=ChainSelector(chain_id))


def renumbering_pMHC(pmhc_paths, output_dir, ref_pdb):

    os.makedirs(output_dir, exist_ok=True)

    for tgt_path in pmhc_paths:

        if not str(tgt_path).lower().endswith('.pdb'):
            if os.path.exists(str(tgt_path)):
                print(f"Skipping non-PDB file: {tgt_path}")
                continue
        else:
            fname = os.path.basename(tgt_path)
        
        
        base, ext = os.path.splitext(fname)
        out_path = os.path.join(output_dir, f"{base}_renum{ext}")
        try:
            parser_pdb = PDBParser(QUIET=True)
            struct = parser_pdb.get_structure('TGT', tgt_path)
            tmp = os.path.join(output_dir, f'.tmp_{fname}')
            best_chain = None
            best_score = -1.0
            best_a1 = best_a2 = None
            for chain in struct[0]:
                extract_chain(tgt_path, tmp, chain.id)
                score, a1, a2 = run_tmalign_with_score(ref_pdb, tmp)

                if score > best_score:
                    best_score, best_chain, best_a1, best_a2 = score, chain.id, a1, a2

            print(f"Best chain for {fname}: {best_chain} with TM-score {best_score:.4f}")
            os.remove(tmp)

            if best_chain is None:
                raise RuntimeError('No suitable chain found for alignment')
            mapping = extract_residue_mapping_with_pdb_numbers(
                ref_pdb, tgt_path, best_a1, best_a2
            )
            mapping = {k:v for k,v in mapping.items() if k[0] == best_chain}
            # fill C-terminus
            parser2 = PDBParser(QUIET=True)
            s2 = parser2.get_structure('TMP2', tgt_path)
            ch = s2[0][best_chain]
            all_res = [r for r in ch if r.id[0] == ' ' and _is_backbone_complete(r)]
            print(_is_backbone_complete(all_res[0]))
            existing = {num for num, cid, _ in mapping.values() if cid == best_chain}
            mx = max(existing) if existing else 0
            for r in all_res:
                key = (best_chain, r.id)
                if key not in mapping:
                    mx += 1
                    mapping[key] = (mx, best_chain, ' ')

            renumber_target_pdb(tgt_path, out_path, mapping)
            print(f"Renumbered chain {best_chain} of {fname} -> {out_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to process {fname}: {e}")
