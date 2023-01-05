import logging
import torch
from ab_rmsd.utils.protein import parsers, constants
from Bio.PDB import Polypeptide


def _aa_tensor_to_sequence(aa):
    return "".join([Polypeptide.index_to_one(a.item()) for a in aa.flatten()])


def _label_heavy_chain_cdr(data, seq_map, max_cdr3_length=30):
    if data is None or seq_map is None:
        return data, seq_map

    # Add CDR labels
    cdr_flag = torch.zeros_like(data["aa"])
    for position, idx in seq_map.items():
        resseq = position[1]
        cdr_type = constants.ChothiaCDRRange.to_cdr("H", resseq)
        if cdr_type is not None:
            cdr_flag[idx] = cdr_type
    data["cdr_flag"] = cdr_flag

    # Add CDR sequence annotations
    data["H1_seq"] = _aa_tensor_to_sequence(data["aa"][cdr_flag == constants.CDR.H1])
    data["H2_seq"] = _aa_tensor_to_sequence(data["aa"][cdr_flag == constants.CDR.H2])
    data["H3_seq"] = _aa_tensor_to_sequence(data["aa"][cdr_flag == constants.CDR.H3])

    cdr3_length = (cdr_flag == constants.CDR.H3).sum().item()
    # Remove too long CDR3
    if cdr3_length > max_cdr3_length:
        cdr_flag[cdr_flag == constants.CDR.H3] = 0
        logging.warning(f"CDR-H3 too long {cdr3_length}. Removed.")
        return None, None

    # Filter: ensure CDR3 exists
    if cdr3_length == 0:
        logging.warning("No CDR-H3 found in the heavy chain.")
        return None, None

    return data, seq_map


def _label_light_chain_cdr(data, seq_map, max_cdr3_length=30):
    if data is None or seq_map is None:
        return data, seq_map
    cdr_flag = torch.zeros_like(data["aa"])
    for position, idx in seq_map.items():
        resseq = position[1]
        cdr_type = constants.ChothiaCDRRange.to_cdr("L", resseq)
        if cdr_type is not None:
            cdr_flag[idx] = cdr_type
    data["cdr_flag"] = cdr_flag

    data["L1_seq"] = _aa_tensor_to_sequence(data["aa"][cdr_flag == constants.CDR.L1])
    data["L2_seq"] = _aa_tensor_to_sequence(data["aa"][cdr_flag == constants.CDR.L2])
    data["L3_seq"] = _aa_tensor_to_sequence(data["aa"][cdr_flag == constants.CDR.L3])

    cdr3_length = (cdr_flag == constants.CDR.L3).sum().item()
    # Remove too long CDR3
    if cdr3_length > max_cdr3_length:
        cdr_flag[cdr_flag == constants.CDR.L3] = 0
        logging.warning(f"CDR-L3 too long {cdr3_length}. Removed.")
        return None, None

    # Ensure CDR3 exists
    if cdr3_length == 0:
        logging.warning("No CDRs found in the light chain.")
        return None, None

    return data, seq_map
