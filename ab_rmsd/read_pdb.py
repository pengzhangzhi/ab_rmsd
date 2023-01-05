import logging
from Bio import PDB
from Bio.PDB import PDBExceptions
from ab_rmsd.utils.label_chain import _label_heavy_chain_cdr, _label_light_chain_cdr
from ab_rmsd.utils.protein import parsers
from ab_rmsd.utils.protein.constants import CDR, CDRID2CDR, NONCDRID
from Bio.PDB import Model


def preprocess_antibody_structure(model:Model, H_id, L_id, id=None):
    """
    process a antibody model

    Args:
        pdb_path (str):
        H_id (str): heavy chain id
        L_id (str): light chain id
        id (str): id of the antibody structure,
    Raises:
        ValueError: if can't parse heavy and light chain from the pdb file

    Returns:
        dict: heavy and light chain information.
    """

    all_chain_ids = [c.id for c in model]

    parsed = {
        "id": id,
        "heavy": None,
        "heavy_seqmap": None,
        "light": None,
        "light_seqmap": None,
        "antigen": None,
        "antigen_seqmap": None,
    }
    try:
        if H_id in all_chain_ids:
            (parsed["heavy"], parsed["heavy_seqmap"]) = _label_heavy_chain_cdr(
                *parsers.parse_biopython_structure(
                    model[H_id], max_resseq=113  # Chothia, end of Heavy chain Fv
                )
            )
            parsed["heavy"]["select"] = {
                CDRID2CDR[flag]: (parsed["heavy"]["cdr_flag"] == flag)
                for flag in [CDR.H1, CDR.H2, CDR.H3]
            }
            parsed["heavy"]["select"].update(
                {"fv-H": (parsed["heavy"]["cdr_flag"] == NONCDRID)}
            )

        if L_id in all_chain_ids:
            (parsed["light"], parsed["light_seqmap"]) = _label_light_chain_cdr(
                *parsers.parse_biopython_structure(
                    model[L_id], max_resseq=106  # Chothia, end of Light chain Fv
                )
            )

            parsed["light"]["select"] = {
                CDRID2CDR[flag]: (parsed["light"]["cdr_flag"] == flag)
                for flag in [CDR.L1, CDR.L2, CDR.L3]
            }
            parsed["light"]["select"].update(
                {"fv-L": (parsed["light"]["cdr_flag"] == NONCDRID)}
            )

        if parsed["heavy"] is None and parsed["light"] is None:
            raise ValueError(
                f"Neither valid antibody H-chain or L-chain is found. "
                f'Please ensure that the chain id of heavy chain is "{H_id}" '
                f'and the id of the light chain is "{L_id}".'
            )

        ag_chain_ids = [cid for cid in all_chain_ids if cid not in (H_id, L_id)]
        if len(ag_chain_ids) > 0:
            chains = [model[c] for c in ag_chain_ids]
            (
                parsed["antigen"],
                parsed["antigen_seqmap"],
            ) = parsers.parse_biopython_structure(chains)

    except (
        PDBExceptions.PDBConstructionException,
        parsers.ParsingException,
        KeyError,
        ValueError,
    ) as e:
        logging.warning("[{}] {}: {}".format(id, e.__class__.__name__, str(e)))
        return None

    return parsed



