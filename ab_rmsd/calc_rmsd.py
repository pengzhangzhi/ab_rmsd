from ab_rmsd.ab_number import renumber
from ab_rmsd.read_pdb import preprocess_antibody_structure
from ab_rmsd.utils.utility import MergeChains, exists, exist_key,PDBParseError, save_pdb
from .superimpose import KabschRMSD
import torch



class AntibodyRMSD:
    def __init__(
        self,
    ):
        self.nb_atoms = 3  # num of atoms selected for superimpose
        self.kabsch_rmsd = KabschRMSD()

    def __call__(self, pred_antibody, native_antibody, superimpose_pred=False):
        """
        calculate the rmsd between two antibodys.
        if superimpose_pred is True, the predicted coordinates will be superimposed to the native coordinates.

        Args:
            pred_antibody (dict): 
            native_antibody (dict): 
            superimpose_pred (bool, optional): superimpose the predicted coordinates. Defaults to True.

        Returns:
            dict: rmsd results of each regions (CDRH1...3, CDRL1...3)
        """
        pred_coord_list, native_coord_list = [], []
        exist_heavy = exist_key(pred_antibody, "heavy") and exist_key(
            native_antibody, "heavy"
        )
        exist_light = exist_key(pred_antibody, "light") and exist_key(
            native_antibody, "light"
        )
        chain_list = []
        if exist_heavy:
            chain_list.append("heavy")
        if exist_light:
            chain_list.append("light")

        for chain in chain_list:
            native_coord = native_antibody[chain]["pos_heavyatom"][
                ..., : self.nb_atoms, :
            ]
            pred_coord = pred_antibody[chain]["pos_heavyatom"][..., : self.nb_atoms, :]
            

            pred_coord_list.append(pred_coord)
            native_coord_list.append(native_coord)

        pred_coord = torch.cat(pred_coord_list, dim=0).reshape(-1, 3)  # [nb_atoms, 3]
        native_coord = torch.cat(native_coord_list, dim=0).reshape(
            -1, 3
        )  # [nb_atoms, 3]
        rot, trans = self.kabsch_rmsd.calc_superimpose_transformation(
            native_coord, pred_coord
        )

        if superimpose_pred:
            for chain in chain_list:
                coord = pred_antibody[chain]["pos_heavyatom"]
                shape = coord.shape
                superimposed = self.kabsch_rmsd.apply_transformation(
                    rot, trans, coord.reshape(-1, 3)
                )
                pred_antibody[chain]["pos_heavyatom"] = superimposed.reshape(*shape)
            
        res_list = []
        for i, chain in enumerate(chain_list):
            pred_coord, native_coord = pred_coord_list[i], native_coord_list[i]
            shape = pred_coord.shape
            superimposed_pred_coord = self.kabsch_rmsd.apply_transformation(
                rot, trans, pred_coord.reshape(-1, 3)
            )
            superimposed_pred_coord = superimposed_pred_coord.reshape(*shape)
            chain_res = self._calc_region_rmsd(
                native_antibody[chain]["select"], superimposed_pred_coord, native_coord
            )
            res_list.append(chain_res)

        res = {}
        for i in range(len(res_list)):
            res.update(res_list[i])
        return res

    def _calc_region_rmsd(self, select_mask_dict, pred_coord, native_coord):
        res_dict = {}
        for region_id, region_mask in select_mask_dict.items():
            pred_region_coord = pred_coord[region_mask].reshape(-1, 3)
            native_region_coord = native_coord[region_mask].reshape(-1, 3)
            rmsd = self.kabsch_rmsd(pred_region_coord, native_region_coord)
            res_dict[region_id] = rmsd
        return res_dict

def calc_ab_rmsd(pred_path, native_path):
    """
    calculate the rmsd between two antibodys.

    Args:
        pred_path (str): 
        native_path (str): 
    """
    pred_ab = parse_pdb(pred_path)
    native_ab = parse_pdb(native_path)
    rmsd = AntibodyRMSD()(pred_ab,native_ab)
    return rmsd

def parse_pdb(pred_path):
    """
    parse a pdb file to protein dict.
    the pdb file must contain at least the heavy chain.
    
    Args:
        pred_path (str): 

    Raises:
        ValueError: No heavy chain found in pdb file.

    Returns:
        dict: dict of heavy and light chain info.
    """
    model, heavy_chains, light_chains, other_chains = renumber(pred_path)
    if len(heavy_chains) == 0: raise PDBParseError("No heavy chain found in pdb file, path: {}".format(pred_path))
    h_id = heavy_chains[0]
    l_id = light_chains[0] if len(light_chains) != 0 else None
    
    ab_dict = preprocess_antibody_structure(
        model, h_id, l_id
    )
    return ab_dict
