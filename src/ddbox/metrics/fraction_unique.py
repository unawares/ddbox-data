from ddbox.metrics.utils import get_molecule_from_smiles_if_valid_or_none
from ddbox.registries import metrics_regstry
from rdkit import Chem


@metrics_regstry.register('fraction_unique')
def fraction_unique(smiles_list, k=None):
    if k is not None:
        smiles_list = smiles_list[:k]
    canonic = set()
    for smiles in smiles_list:
        molecule = get_molecule_from_smiles_if_valid_or_none(smiles)
        if molecule is not None:
            canonic.add(Chem.MolToSmiles(molecule))
    if None in canonic:
        raise Exception("canonic set contains an invalid smiles")
    return len(canonic) / len(smiles_list)
