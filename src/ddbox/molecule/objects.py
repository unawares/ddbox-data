from numbers import Number
from typing import *

from ddbox.molecule.descriptors import compute_descriptors
from ddbox.utils import json_to_pretty_str
from rdkit import Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

TMoleculeDescriptor = TypeVar("TMoleculeDescriptor", bound="MoleculeDescriptor")


class MoleculeDescriptor:

    # https://datagrok.ai/help/domains/chem/descriptors

    BalabanJ: Number = None
    BertzCT: Number = None
    Chi0: Number = None
    Chi0n: Number = None
    Chi0v: Number = None
    Chi1: Number = None
    Chi1n: Number = None
    Chi1v: Number = None
    Chi2n: Number = None
    Chi2v: Number = None
    Chi3n: Number = None
    Chi3v: Number = None
    Chi4n: Number = None
    Chi4v: Number = None
    EState_VSA1: Number = None
    EState_VSA10: Number = None
    EState_VSA11: Number = None
    EState_VSA2: Number = None
    EState_VSA3: Number = None
    EState_VSA4: Number = None
    EState_VSA5: Number = None
    EState_VSA6: Number = None
    EState_VSA7: Number = None
    EState_VSA8: Number = None
    EState_VSA9: Number = None
    ExactMolWt: Number = None
    FpDensityMorgan1: Number = None
    FpDensityMorgan2: Number = None
    FpDensityMorgan3: Number = None
    FractionCSP3: Number = None
    HallKierAlpha: Number = None
    HeavyAtomCount: Number = None
    HeavyAtomMolWt: Number = None
    Ipc: Number = None
    Kappa1: Number = None
    Kappa2: Number = None
    Kappa3: Number = None
    LabuteASA: Number = None
    MaxAbsEStateIndex: Number = None
    MaxAbsPartialCharge: Number = None
    MaxEStateIndex: Number = None
    MaxPartialCharge: Number = None
    MinAbsEStateIndex: Number = None
    MinAbsPartialCharge: Number = None
    MinEStateIndex: Number = None
    MinPartialCharge: Number = None
    MolLogP: Number = None
    MolMR: Number = None
    MolWt: Number = None
    NHOHCount: Number = None
    NOCount: Number = None
    NumAliphaticCarbocycles: Number = None
    NumAliphaticHeterocycles: Number = None
    NumAliphaticRings: Number = None
    NumAromaticCarbocycles: Number = None
    NumAromaticHeterocycles: Number = None
    NumAromaticRings: Number = None
    NumHAcceptors: Number = None
    NumHDonors: Number = None
    NumHeteroatoms: Number = None
    NumRadicalElectrons: Number = None
    NumRotatableBonds: Number = None
    NumSaturatedCarbocycles: Number = None
    NumSaturatedHeterocycles: Number = None
    NumSaturatedRings: Number = None
    NumValenceElectrons: Number = None
    PEOE_VSA1: Number = None
    PEOE_VSA10: Number = None
    PEOE_VSA11: Number = None
    PEOE_VSA12: Number = None
    PEOE_VSA13: Number = None
    PEOE_VSA14: Number = None
    PEOE_VSA2: Number = None
    PEOE_VSA3: Number = None
    PEOE_VSA4: Number = None
    PEOE_VSA5: Number = None
    PEOE_VSA6: Number = None
    PEOE_VSA7: Number = None
    PEOE_VSA8: Number = None
    PEOE_VSA9: Number = None
    RingCount: Number = None
    SMR_VSA1: Number = None
    SMR_VSA10: Number = None
    SMR_VSA2: Number = None
    SMR_VSA3: Number = None
    SMR_VSA4: Number = None
    SMR_VSA5: Number = None
    SMR_VSA6: Number = None
    SMR_VSA7: Number = None
    SMR_VSA8: Number = None
    SMR_VSA9: Number = None
    SlogP_VSA1: Number = None
    SlogP_VSA10: Number = None
    SlogP_VSA11: Number = None
    SlogP_VSA12: Number = None
    SlogP_VSA2: Number = None
    SlogP_VSA3: Number = None
    SlogP_VSA4: Number = None
    SlogP_VSA5: Number = None
    SlogP_VSA6: Number = None
    SlogP_VSA7: Number = None
    SlogP_VSA8: Number = None
    SlogP_VSA9: Number = None
    TPSA: Number = None
    VSA_EState1: Number = None
    VSA_EState10: Number = None
    VSA_EState2: Number = None
    VSA_EState3: Number = None
    VSA_EState4: Number = None
    VSA_EState5: Number = None
    VSA_EState6: Number = None
    VSA_EState7: Number = None
    VSA_EState8: Number = None
    VSA_EState9: Number = None
    fr_Al_COO: Number = None
    fr_Al_OH: Number = None
    fr_Al_OH_noTert: Number = None
    fr_ArN: Number = None
    fr_Ar_COO: Number = None
    fr_Ar_N: Number = None
    fr_Ar_NH: Number = None
    fr_Ar_OH: Number = None
    fr_COO: Number = None
    fr_COO2: Number = None
    fr_C_O: Number = None
    fr_C_O_noCOO: Number = None
    fr_C_S: Number = None
    fr_HOCCN: Number = None
    fr_Imine: Number = None
    fr_NH0: Number = None
    fr_NH1: Number = None
    fr_NH2: Number = None
    fr_N_O: Number = None
    fr_Ndealkylation1: Number = None
    fr_Ndealkylation2: Number = None
    fr_Nhpyrrole: Number = None
    fr_SH: Number = None
    fr_aldehyde: Number = None
    fr_alkyl_carbamate: Number = None
    fr_alkyl_halide: Number = None
    fr_allylic_oxid: Number = None
    fr_amide: Number = None
    fr_amidine: Number = None
    fr_aniline: Number = None
    fr_aryl_methyl: Number = None
    fr_azide: Number = None
    fr_azo: Number = None
    fr_barbitur: Number = None
    fr_benzene: Number = None
    fr_benzodiazepine: Number = None
    fr_bicyclic: Number = None
    fr_diazo: Number = None
    fr_dihydropyridine: Number = None
    fr_epoxide: Number = None
    fr_ester: Number = None
    fr_ether: Number = None
    fr_furan: Number = None
    fr_guanido: Number = None
    fr_halogen: Number = None
    fr_hdrzine: Number = None
    fr_hdrzone: Number = None
    fr_imidazole: Number = None
    fr_imide: Number = None
    fr_isocyan: Number = None
    fr_isothiocyan: Number = None
    fr_ketone: Number = None
    fr_ketone_Topliss: Number = None
    fr_lactam: Number = None
    fr_lactone: Number = None
    fr_methoxy: Number = None
    fr_morpholine: Number = None
    fr_nitrile: Number = None
    fr_nitro: Number = None
    fr_nitro_arom: Number = None
    fr_nitro_arom_nonortho: Number = None
    fr_nitroso: Number = None
    fr_oxazole: Number = None
    fr_oxime: Number = None
    fr_para_hydroxylation: Number = None
    fr_phenol: Number = None
    fr_phenol_noOrthoHbond: Number = None
    fr_phos_acid: Number = None
    fr_phos_ester: Number = None
    fr_piperdine: Number = None
    fr_piperzine: Number = None
    fr_priamide: Number = None
    fr_prisulfonamd: Number = None
    fr_pyridine: Number = None
    fr_quatN: Number = None
    fr_sulfide: Number = None
    fr_sulfonamd: Number = None
    fr_sulfone: Number = None
    fr_term_acetylene: Number = None
    fr_tetrazole: Number = None
    fr_thiazole: Number = None
    fr_thiocyan: Number = None
    fr_thiophene: Number = None
    fr_unbrch_alkane: Number = None
    fr_urea: Number = None
    qed: Number = None

    def __init__(self) -> None:
        pass

    @property
    def _descriptors(self) -> List[str]:
        descriptors = [v for v in self.__class__.__dict__.keys() if not v.startswith('__') and not v.startswith('_')]
        return descriptors

    def _to_json_data(self) -> Mapping[str, Any]:
        return {descriptor: self.__getattribute__(descriptor) for descriptor in self._descriptors}


TMolecule = TypeVar("TMolecule", bound="Molecule")


class Molecule:

    _inchi_key: str = None
    _inchi: str = None
    _smiles: str = None

    mol: Chem.rdchem.Mol
    desciptor: MoleculeDescriptor
    is_fetched: bool = False

    def __init__(self) -> None:
        self.desciptor = MoleculeDescriptor()

    @property
    def inchi_key(self) -> str:
        if self.is_valid:
            return Chem.MolToInchiKey(self.mol)
        return self._inchi_key

    @property
    def inchi(self) -> str:
        if self.is_valid:
            return Chem.MolToInchi(self.mol)
        return self._inchi

    @property
    def smiles(self) -> str:
        if self.is_valid:
            return Chem.MolToSmiles(self.mol)
        return self._smiles

    @property
    def is_valid(self) -> str:
        if self.mol is not None:
            try:
                Chem.SanitizeMol(self.mol)
                return True
            except ValueError:
                pass
        return False

    def _fetch_descriptors(self) -> TMolecule:
        if self.is_valid:
            mol_descriptor_calculator = MolecularDescriptorCalculator(self.desciptor._descriptors)
            values = list(mol_descriptor_calculator.CalcDescriptors(self.mol))
            for index in range(len(self.desciptor._descriptors)):
                self.desciptor.__setattr__(self.desciptor._descriptors[index], values[index])
        return self

    def _override_data(self):
        self._inchi_key = self.inchi_key
        self._inchi = self.inchi
        self._smiles = self.smiles
        self._fetch_descriptors()

    def from_inchi(self, inchi) -> TMolecule:
        self._inchi = inchi
        self.mol = Chem.MolFromInchi(inchi)
        self._override_data()
        return self

    def from_smiles(self, smiles) -> TMolecule:
        self._smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        self._override_data()
        return self

    def __str__(self) -> str:
        data = {
            'inchi_key': self.inchi_key,
            'inchi': self.inchi,
            'smiles': self.smiles,
            'descriptors': self.desciptor._to_json_data()
        }
        return json_to_pretty_str(data)
