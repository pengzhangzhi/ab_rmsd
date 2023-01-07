import torch
import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from .protein.constants import AA, Fragment, restype_to_heavyatom_names


def exists(x):
    return x is not None


def exist_key(d, k):
    return d.get(k) is not None

class PDBParseError(Exception):
    pass

def save_pdb(data, path=None):
    """
    Args:
        data:   A dict that contains: `chain_nb`, `chain_id`, `aa`, `resseq`, `icode`,
                `pos_heavyatom`, `mask_heavyatom`.
    """

    def _mask_select(v, mask):
        if isinstance(v, str):
            return ''.join([s for i, s in enumerate(v) if mask[i]])
        elif isinstance(v, list):
            return [s for i, s in enumerate(v) if mask[i]]
        elif isinstance(v, torch.Tensor):
            return v[mask]
        else:
            return v

    def _build_chain(builder, aa_ch, pos_heavyatom_ch, mask_heavyatom_ch, chain_id_ch, resseq_ch, icode_ch):
        builder.init_chain(chain_id_ch[0])
        builder.init_seg('    ')

        for aa_res, pos_allatom_res, mask_allatom_res, resseq_res, icode_res in \
            zip(aa_ch, pos_heavyatom_ch, mask_heavyatom_ch, resseq_ch, icode_ch):
            if not AA.is_aa(aa_res.item()): 
                print('[Warning] Unknown amino acid type at %d%s: %r' % (resseq_res.item(), icode_res, aa_res.item()))
                continue
            restype = AA(aa_res.item())
            builder.init_residue(
                resname = str(restype),
                field = ' ',
                resseq = resseq_res.item(),
                icode = icode_res,
            )

            for i, atom_name in enumerate(restype_to_heavyatom_names[restype]):
                if atom_name == '': continue    # No expected atom
                if (~mask_allatom_res[i]).any(): continue     # Atom is missing
                if len(atom_name) == 1: fullname = ' %s  ' % atom_name
                elif len(atom_name) == 2: fullname = ' %s ' % atom_name
                elif len(atom_name) == 3: fullname = ' %s' % atom_name
                else: fullname = atom_name # len == 4
                builder.init_atom(atom_name, pos_allatom_res[i].tolist(), 0.0, 1.0, ' ', fullname,)

    warnings.simplefilter('ignore', BiopythonWarning)
    builder = StructureBuilder()
    builder.init_structure(0)
    builder.init_model(0)

    unique_chain_nb = data['chain_nb'].unique().tolist()
    for ch_nb in unique_chain_nb:
        mask = (data['chain_nb'] == ch_nb)
        aa = _mask_select(data['aa'], mask)
        pos_heavyatom = _mask_select(data['pos_heavyatom'], mask)
        mask_heavyatom = _mask_select(data['mask_heavyatom'], mask)
        chain_id = _mask_select(data['chain_id'], mask)
        resseq = _mask_select(data['resseq'], mask)
        icode = _mask_select(data['icode'], mask)

        _build_chain(builder, aa, pos_heavyatom, mask_heavyatom, chain_id, resseq, icode)
    
    structure = builder.get_structure()
    if path is not None:
        io = PDBIO()
        io.set_structure(structure)
        io.save(path)
    return structure



class MergeChains(object):

    def __init__(self):
        super().__init__()

    def assign_chain_number_(self, data_list):
        chains = set()
        for data in data_list:
            chains.update(data['chain_id'])
        # assign a chain id to a number (0, 1, 2, 3, ...)
        chains = {c: i for i, c in enumerate(chains)}

        for data in data_list:
            data['chain_nb'] = torch.LongTensor([
                chains[c] for c in data['chain_id']
            ])

    def _data_attr(self, data, name):
        if name in ('generate_flag', 'anchor_flag') and name not in data:
            return torch.zeros(data['aa'].shape, dtype=torch.bool)
        else:
            return data[name]

    def __call__(self, structure):
        """ merge antigen, heavy and light chains data into one tensor."""
        data_list = []
        if structure['heavy'] is not None:
            structure['heavy']['fragment_type'] = torch.full_like(
                structure['heavy']['aa'],
                fill_value = Fragment.Heavy,
            )
            data_list.append(structure['heavy'])

        if structure['light'] is not None:
            structure['light']['fragment_type'] = torch.full_like(
                structure['light']['aa'],
                fill_value = Fragment.Light,
            )
            data_list.append(structure['light'])

        if structure['antigen'] is not None:
            structure['antigen']['fragment_type'] = torch.full_like(
                structure['antigen']['aa'],
                fill_value = Fragment.Antigen,
            )
            structure['antigen']['cdr_flag'] = torch.zeros_like(
                structure['antigen']['aa'],
            )
            data_list.append(structure['antigen'])

        self.assign_chain_number_(data_list)

        list_props = {
            'chain_id': [],
            'icode': [],
        }
        tensor_props = {
            'chain_nb': [],
            'resseq': [],
            'res_nb': [],
            'aa': [],
            'pos_heavyatom': [],
            'mask_heavyatom': [],
            'generate_flag': [],
            'cdr_flag': [],
            'anchor_flag': [],
            'fragment_type': [],
        }

        for data in data_list:
            for k in list_props.keys():
                list_props[k].append(self._data_attr(data, k))
            for k in tensor_props.keys():
                tensor_props[k].append(self._data_attr(data, k))

        list_props = {k: sum(v, start=[]) for k, v in list_props.items()}
        tensor_props = {k: torch.cat(v, dim=0) for k, v in tensor_props.items()}
        data_out = {
            **list_props,
            **tensor_props,
        }
        return data_out

