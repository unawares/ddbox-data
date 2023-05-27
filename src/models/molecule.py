from datetime import date, datetime
from typing import *

from app.configs import settings
from models.base import Base
from sqlalchemy import URL, BigInteger, Column, Date, DateTime, Float, ForeignKey, Integer, String, Table
from sqlalchemy.orm import relationship

molecule_tag_m2m = Table('molecule_tag_m2m', Base.metadata,
    Column('molecule_id', settings.SQLALCHEMY_ID_TYPE, ForeignKey('molecules.id'), index=True),
    Column('tag_id', settings.SQLALCHEMY_ID_TYPE, ForeignKey('tags.id'), index=True)
 )


TTagModel = TypeVar('TTagModel', bound='TagModel')


class TagModel(Base):
    __tablename__ = 'tags'

    id = Column(settings.SQLALCHEMY_ID_TYPE, primary_key=True, index=True)
    name = Column(String, unique=True)

    def __init__(
            self: TTagModel,
            name: str,
    ):
        self.name = name

    def __repr__(self: TTagModel):
        return "<TagModel(id='%s', name='%s')>" % (
            self.id,
            self.name,
        )


TMoleculeModel = TypeVar('TMoleculeModel', bound='MoleculeModel')


class MoleculeModel(Base):
    __tablename__ = 'molecules'

    id = Column(settings.SQLALCHEMY_ID_TYPE, primary_key=True, index=True)

    inchi_key = Column(String)
    inchi = Column(String)
    smiles = Column(String)

    BalabanJ = Column(Float)
    BertzCT = Column(Float)
    Chi0 = Column(Float)
    Chi0n = Column(Float)
    Chi0v = Column(Float)
    Chi1 = Column(Float)
    Chi1n = Column(Float)
    Chi1v = Column(Float)
    Chi2n = Column(Float)
    Chi2v = Column(Float)
    Chi3n = Column(Float)
    Chi3v = Column(Float)
    Chi4n = Column(Float)
    Chi4v = Column(Float)
    EState_VSA1 = Column(Float)
    EState_VSA10 = Column(Float)
    EState_VSA11 = Column(Float)
    EState_VSA2 = Column(Float)
    EState_VSA3 = Column(Float)
    EState_VSA4 = Column(Float)
    EState_VSA5 = Column(Float)
    EState_VSA6 = Column(Float)
    EState_VSA7 = Column(Float)
    EState_VSA8 = Column(Float)
    EState_VSA9 = Column(Float)
    ExactMolWt = Column(Float)
    FpDensityMorgan1 = Column(Float)
    FpDensityMorgan2 = Column(Float)
    FpDensityMorgan3 = Column(Float)
    FractionCSP3 = Column(Float)
    HallKierAlpha = Column(Float)
    HeavyAtomCount = Column(Float)
    HeavyAtomMolWt = Column(Float)
    Ipc = Column(Float)
    Kappa1 = Column(Float)
    Kappa2 = Column(Float)
    Kappa3 = Column(Float)
    LabuteASA = Column(Float)
    MaxAbsEStateIndex = Column(Float)
    MaxAbsPartialCharge = Column(Float)
    MaxEStateIndex = Column(Float)
    MaxPartialCharge = Column(Float)
    MinAbsEStateIndex = Column(Float)
    MinAbsPartialCharge = Column(Float)
    MinEStateIndex = Column(Float)
    MinPartialCharge = Column(Float)
    MolLogP = Column(Float)
    MolMR = Column(Float)
    MolWt = Column(Float)
    NHOHCount = Column(Float)
    NOCount = Column(Float)
    NumAliphaticCarbocycles = Column(Float)
    NumAliphaticHeterocycles = Column(Float)
    NumAliphaticRings = Column(Float)
    NumAromaticCarbocycles = Column(Float)
    NumAromaticHeterocycles = Column(Float)
    NumAromaticRings = Column(Float)
    NumHAcceptors = Column(Float)
    NumHDonors = Column(Float)
    NumHeteroatoms = Column(Float)
    NumRadicalElectrons = Column(Float)
    NumRotatableBonds = Column(Float)
    NumSaturatedCarbocycles = Column(Float)
    NumSaturatedHeterocycles = Column(Float)
    NumSaturatedRings = Column(Float)
    NumValenceElectrons = Column(Float)
    PEOE_VSA1 = Column(Float)
    PEOE_VSA10 = Column(Float)
    PEOE_VSA11 = Column(Float)
    PEOE_VSA12 = Column(Float)
    PEOE_VSA13 = Column(Float)
    PEOE_VSA14 = Column(Float)
    PEOE_VSA2 = Column(Float)
    PEOE_VSA3 = Column(Float)
    PEOE_VSA4 = Column(Float)
    PEOE_VSA5 = Column(Float)
    PEOE_VSA6 = Column(Float)
    PEOE_VSA7 = Column(Float)
    PEOE_VSA8 = Column(Float)
    PEOE_VSA9 = Column(Float)
    RingCount = Column(Float)
    SMR_VSA1 = Column(Float)
    SMR_VSA10 = Column(Float)
    SMR_VSA2 = Column(Float)
    SMR_VSA3 = Column(Float)
    SMR_VSA4 = Column(Float)
    SMR_VSA5 = Column(Float)
    SMR_VSA6 = Column(Float)
    SMR_VSA7 = Column(Float)
    SMR_VSA8 = Column(Float)
    SMR_VSA9 = Column(Float)
    SlogP_VSA1 = Column(Float)
    SlogP_VSA10 = Column(Float)
    SlogP_VSA11 = Column(Float)
    SlogP_VSA12 = Column(Float)
    SlogP_VSA2 = Column(Float)
    SlogP_VSA3 = Column(Float)
    SlogP_VSA4 = Column(Float)
    SlogP_VSA5 = Column(Float)
    SlogP_VSA6 = Column(Float)
    SlogP_VSA7 = Column(Float)
    SlogP_VSA8 = Column(Float)
    SlogP_VSA9 = Column(Float)
    TPSA = Column(Float)
    VSA_EState1 = Column(Float)
    VSA_EState10 = Column(Float)
    VSA_EState2 = Column(Float)
    VSA_EState3 = Column(Float)
    VSA_EState4 = Column(Float)
    VSA_EState5 = Column(Float)
    VSA_EState6 = Column(Float)
    VSA_EState7 = Column(Float)
    VSA_EState8 = Column(Float)
    VSA_EState9 = Column(Float)
    fr_Al_COO = Column(Float)
    fr_Al_OH = Column(Float)
    fr_Al_OH_noTert = Column(Float)
    fr_ArN = Column(Float)
    fr_Ar_COO = Column(Float)
    fr_Ar_N = Column(Float)
    fr_Ar_NH = Column(Float)
    fr_Ar_OH = Column(Float)
    fr_COO = Column(Float)
    fr_COO2 = Column(Float)
    fr_C_O = Column(Float)
    fr_C_O_noCOO = Column(Float)
    fr_C_S = Column(Float)
    fr_HOCCN = Column(Float)
    fr_Imine = Column(Float)
    fr_NH0 = Column(Float)
    fr_NH1 = Column(Float)
    fr_NH2 = Column(Float)
    fr_N_O = Column(Float)
    fr_Ndealkylation1 = Column(Float)
    fr_Ndealkylation2 = Column(Float)
    fr_Nhpyrrole = Column(Float)
    fr_SH = Column(Float)
    fr_aldehyde = Column(Float)
    fr_alkyl_carbamate = Column(Float)
    fr_alkyl_halide = Column(Float)
    fr_allylic_oxid = Column(Float)
    fr_amide = Column(Float)
    fr_amidine = Column(Float)
    fr_aniline = Column(Float)
    fr_aryl_methyl = Column(Float)
    fr_azide = Column(Float)
    fr_azo = Column(Float)
    fr_barbitur = Column(Float)
    fr_benzene = Column(Float)
    fr_benzodiazepine = Column(Float)
    fr_bicyclic = Column(Float)
    fr_diazo = Column(Float)
    fr_dihydropyridine = Column(Float)
    fr_epoxide = Column(Float)
    fr_ester = Column(Float)
    fr_ether = Column(Float)
    fr_furan = Column(Float)
    fr_guanido = Column(Float)
    fr_halogen = Column(Float)
    fr_hdrzine = Column(Float)
    fr_hdrzone = Column(Float)
    fr_imidazole = Column(Float)
    fr_imide = Column(Float)
    fr_isocyan = Column(Float)
    fr_isothiocyan = Column(Float)
    fr_ketone = Column(Float)
    fr_ketone_Topliss = Column(Float)
    fr_lactam = Column(Float)
    fr_lactone = Column(Float)
    fr_methoxy = Column(Float)
    fr_morpholine = Column(Float)
    fr_nitrile = Column(Float)
    fr_nitro = Column(Float)
    fr_nitro_arom = Column(Float)
    fr_nitro_arom_nonortho = Column(Float)
    fr_nitroso = Column(Float)
    fr_oxazole = Column(Float)
    fr_oxime = Column(Float)
    fr_para_hydroxylation = Column(Float)
    fr_phenol = Column(Float)
    fr_phenol_noOrthoHbond = Column(Float)
    fr_phos_acid = Column(Float)
    fr_phos_ester = Column(Float)
    fr_piperdine = Column(Float)
    fr_piperzine = Column(Float)
    fr_priamide = Column(Float)
    fr_prisulfonamd = Column(Float)
    fr_pyridine = Column(Float)
    fr_quatN = Column(Float)
    fr_sulfide = Column(Float)
    fr_sulfonamd = Column(Float)
    fr_sulfone = Column(Float)
    fr_term_acetylene = Column(Float)
    fr_tetrazole = Column(Float)
    fr_thiazole = Column(Float)
    fr_thiocyan = Column(Float)
    fr_thiophene = Column(Float)
    fr_unbrch_alkane = Column(Float)
    fr_urea = Column(Float)
    qed = Column(Float)

    tags = relationship('TagModel', secondary=molecule_tag_m2m, backref='molecules')

    def __init__(
            self: TMoleculeModel,
            inchi_key: str,
            inchi: str,
            smiles: str,
            BalabanJ: str,
            BertzCT: str,
            Chi0: str,
            Chi0n: str,
            Chi0v: str,
            Chi1: str,
            Chi1n: str,
            Chi1v: str,
            Chi2n: str,
            Chi2v: str,
            Chi3n: str,
            Chi3v: str,
            Chi4n: str,
            Chi4v: str,
            EState_VSA1: str,
            EState_VSA10: str,
            EState_VSA11: str,
            EState_VSA2: str,
            EState_VSA3: str,
            EState_VSA4: str,
            EState_VSA5: str,
            EState_VSA6: str,
            EState_VSA7: str,
            EState_VSA8: str,
            EState_VSA9: str,
            ExactMolWt: str,
            FpDensityMorgan1: str,
            FpDensityMorgan2: str,
            FpDensityMorgan3: str,
            FractionCSP3: str,
            HallKierAlpha: str,
            HeavyAtomCount: str,
            HeavyAtomMolWt: str,
            Ipc: str,
            Kappa1: str,
            Kappa2: str,
            Kappa3: str,
            LabuteASA: str,
            MaxAbsEStateIndex: str,
            MaxAbsPartialCharge: str,
            MaxEStateIndex: str,
            MaxPartialCharge: str,
            MinAbsEStateIndex: str,
            MinAbsPartialCharge: str,
            MinEStateIndex: str,
            MinPartialCharge: str,
            MolLogP: str,
            MolMR: str,
            MolWt: str,
            NHOHCount: str,
            NOCount: str,
            NumAliphaticCarbocycles: str,
            NumAliphaticHeterocycles: str,
            NumAliphaticRings: str,
            NumAromaticCarbocycles: str,
            NumAromaticHeterocycles: str,
            NumAromaticRings: str,
            NumHAcceptors: str,
            NumHDonors: str,
            NumHeteroatoms: str,
            NumRadicalElectrons: str,
            NumRotatableBonds: str,
            NumSaturatedCarbocycles: str,
            NumSaturatedHeterocycles: str,
            NumSaturatedRings: str,
            NumValenceElectrons: str,
            PEOE_VSA1: str,
            PEOE_VSA10: str,
            PEOE_VSA11: str,
            PEOE_VSA12: str,
            PEOE_VSA13: str,
            PEOE_VSA14: str,
            PEOE_VSA2: str,
            PEOE_VSA3: str,
            PEOE_VSA4: str,
            PEOE_VSA5: str,
            PEOE_VSA6: str,
            PEOE_VSA7: str,
            PEOE_VSA8: str,
            PEOE_VSA9: str,
            RingCount: str,
            SMR_VSA1: str,
            SMR_VSA10: str,
            SMR_VSA2: str,
            SMR_VSA3: str,
            SMR_VSA4: str,
            SMR_VSA5: str,
            SMR_VSA6: str,
            SMR_VSA7: str,
            SMR_VSA8: str,
            SMR_VSA9: str,
            SlogP_VSA1: str,
            SlogP_VSA10: str,
            SlogP_VSA11: str,
            SlogP_VSA12: str,
            SlogP_VSA2: str,
            SlogP_VSA3: str,
            SlogP_VSA4: str,
            SlogP_VSA5: str,
            SlogP_VSA6: str,
            SlogP_VSA7: str,
            SlogP_VSA8: str,
            SlogP_VSA9: str,
            TPSA: str,
            VSA_EState1: str,
            VSA_EState10: str,
            VSA_EState2: str,
            VSA_EState3: str,
            VSA_EState4: str,
            VSA_EState5: str,
            VSA_EState6: str,
            VSA_EState7: str,
            VSA_EState8: str,
            VSA_EState9: str,
            fr_Al_COO: str,
            fr_Al_OH: str,
            fr_Al_OH_noTert: str,
            fr_ArN: str,
            fr_Ar_COO: str,
            fr_Ar_N: str,
            fr_Ar_NH: str,
            fr_Ar_OH: str,
            fr_COO: str,
            fr_COO2: str,
            fr_C_O: str,
            fr_C_O_noCOO: str,
            fr_C_S: str,
            fr_HOCCN: str,
            fr_Imine: str,
            fr_NH0: str,
            fr_NH1: str,
            fr_NH2: str,
            fr_N_O: str,
            fr_Ndealkylation1: str,
            fr_Ndealkylation2: str,
            fr_Nhpyrrole: str,
            fr_SH: str,
            fr_aldehyde: str,
            fr_alkyl_carbamate: str,
            fr_alkyl_halide: str,
            fr_allylic_oxid: str,
            fr_amide: str,
            fr_amidine: str,
            fr_aniline: str,
            fr_aryl_methyl: str,
            fr_azide: str,
            fr_azo: str,
            fr_barbitur: str,
            fr_benzene: str,
            fr_benzodiazepine: str,
            fr_bicyclic: str,
            fr_diazo: str,
            fr_dihydropyridine: str,
            fr_epoxide: str,
            fr_ester: str,
            fr_ether: str,
            fr_furan: str,
            fr_guanido: str,
            fr_halogen: str,
            fr_hdrzine: str,
            fr_hdrzone: str,
            fr_imidazole: str,
            fr_imide: str,
            fr_isocyan: str,
            fr_isothiocyan: str,
            fr_ketone: str,
            fr_ketone_Topliss: str,
            fr_lactam: str,
            fr_lactone: str,
            fr_methoxy: str,
            fr_morpholine: str,
            fr_nitrile: str,
            fr_nitro: str,
            fr_nitro_arom: str,
            fr_nitro_arom_nonortho: str,
            fr_nitroso: str,
            fr_oxazole: str,
            fr_oxime: str,
            fr_para_hydroxylation: str,
            fr_phenol: str,
            fr_phenol_noOrthoHbond: str,
            fr_phos_acid: str,
            fr_phos_ester: str,
            fr_piperdine: str,
            fr_piperzine: str,
            fr_priamide: str,
            fr_prisulfonamd: str,
            fr_pyridine: str,
            fr_quatN: str,
            fr_sulfide: str,
            fr_sulfonamd: str,
            fr_sulfone: str,
            fr_term_acetylene: str,
            fr_tetrazole: str,
            fr_thiazole: str,
            fr_thiocyan: str,
            fr_thiophene: str,
            fr_unbrch_alkane: str,
            fr_urea: str,
            qed: str,
            tags: List[TagModel],
    ):
        self.inchi_key = inchi_key
        self.inchi = inchi
        self.smiles = smiles
        self.BalabanJ = BalabanJ
        self.BertzCT = BertzCT
        self.Chi0 = Chi0
        self.Chi0n = Chi0n
        self.Chi0v = Chi0v
        self.Chi1 = Chi1
        self.Chi1n = Chi1n
        self.Chi1v = Chi1v
        self.Chi2n = Chi2n
        self.Chi2v = Chi2v
        self.Chi3n = Chi3n
        self.Chi3v = Chi3v
        self.Chi4n = Chi4n
        self.Chi4v = Chi4v
        self.EState_VSA1 = EState_VSA1
        self.EState_VSA10 = EState_VSA10
        self.EState_VSA11 = EState_VSA11
        self.EState_VSA2 = EState_VSA2
        self.EState_VSA3 = EState_VSA3
        self.EState_VSA4 = EState_VSA4
        self.EState_VSA5 = EState_VSA5
        self.EState_VSA6 = EState_VSA6
        self.EState_VSA7 = EState_VSA7
        self.EState_VSA8 = EState_VSA8
        self.EState_VSA9 = EState_VSA9
        self.ExactMolWt = ExactMolWt
        self.FpDensityMorgan1 = FpDensityMorgan1
        self.FpDensityMorgan2 = FpDensityMorgan2
        self.FpDensityMorgan3 = FpDensityMorgan3
        self.FractionCSP3 = FractionCSP3
        self.HallKierAlpha = HallKierAlpha
        self.HeavyAtomCount = HeavyAtomCount
        self.HeavyAtomMolWt = HeavyAtomMolWt
        self.Ipc = Ipc
        self.Kappa1 = Kappa1
        self.Kappa2 = Kappa2
        self.Kappa3 = Kappa3
        self.LabuteASA = LabuteASA
        self.MaxAbsEStateIndex = MaxAbsEStateIndex
        self.MaxAbsPartialCharge = MaxAbsPartialCharge
        self.MaxEStateIndex = MaxEStateIndex
        self.MaxPartialCharge = MaxPartialCharge
        self.MinAbsEStateIndex = MinAbsEStateIndex
        self.MinAbsPartialCharge = MinAbsPartialCharge
        self.MinEStateIndex = MinEStateIndex
        self.MinPartialCharge = MinPartialCharge
        self.MolLogP = MolLogP
        self.MolMR = MolMR
        self.MolWt = MolWt
        self.NHOHCount = NHOHCount
        self.NOCount = NOCount
        self.NumAliphaticCarbocycles = NumAliphaticCarbocycles
        self.NumAliphaticHeterocycles = NumAliphaticHeterocycles
        self.NumAliphaticRings = NumAliphaticRings
        self.NumAromaticCarbocycles = NumAromaticCarbocycles
        self.NumAromaticHeterocycles = NumAromaticHeterocycles
        self.NumAromaticRings = NumAromaticRings
        self.NumHAcceptors = NumHAcceptors
        self.NumHDonors = NumHDonors
        self.NumHeteroatoms = NumHeteroatoms
        self.NumRadicalElectrons = NumRadicalElectrons
        self.NumRotatableBonds = NumRotatableBonds
        self.NumSaturatedCarbocycles = NumSaturatedCarbocycles
        self.NumSaturatedHeterocycles = NumSaturatedHeterocycles
        self.NumSaturatedRings = NumSaturatedRings
        self.NumValenceElectrons = NumValenceElectrons
        self.PEOE_VSA1 = PEOE_VSA1
        self.PEOE_VSA10 = PEOE_VSA10
        self.PEOE_VSA11 = PEOE_VSA11
        self.PEOE_VSA12 = PEOE_VSA12
        self.PEOE_VSA13 = PEOE_VSA13
        self.PEOE_VSA14 = PEOE_VSA14
        self.PEOE_VSA2 = PEOE_VSA2
        self.PEOE_VSA3 = PEOE_VSA3
        self.PEOE_VSA4 = PEOE_VSA4
        self.PEOE_VSA5 = PEOE_VSA5
        self.PEOE_VSA6 = PEOE_VSA6
        self.PEOE_VSA7 = PEOE_VSA7
        self.PEOE_VSA8 = PEOE_VSA8
        self.PEOE_VSA9 = PEOE_VSA9
        self.RingCount = RingCount
        self.SMR_VSA1 = SMR_VSA1
        self.SMR_VSA10 = SMR_VSA10
        self.SMR_VSA2 = SMR_VSA2
        self.SMR_VSA3 = SMR_VSA3
        self.SMR_VSA4 = SMR_VSA4
        self.SMR_VSA5 = SMR_VSA5
        self.SMR_VSA6 = SMR_VSA6
        self.SMR_VSA7 = SMR_VSA7
        self.SMR_VSA8 = SMR_VSA8
        self.SMR_VSA9 = SMR_VSA9
        self.SlogP_VSA1 = SlogP_VSA1
        self.SlogP_VSA10 = SlogP_VSA10
        self.SlogP_VSA11 = SlogP_VSA11
        self.SlogP_VSA12 = SlogP_VSA12
        self.SlogP_VSA2 = SlogP_VSA2
        self.SlogP_VSA3 = SlogP_VSA3
        self.SlogP_VSA4 = SlogP_VSA4
        self.SlogP_VSA5 = SlogP_VSA5
        self.SlogP_VSA6 = SlogP_VSA6
        self.SlogP_VSA7 = SlogP_VSA7
        self.SlogP_VSA8 = SlogP_VSA8
        self.SlogP_VSA9 = SlogP_VSA9
        self.TPSA = TPSA
        self.VSA_EState1 = VSA_EState1
        self.VSA_EState10 = VSA_EState10
        self.VSA_EState2 = VSA_EState2
        self.VSA_EState3 = VSA_EState3
        self.VSA_EState4 = VSA_EState4
        self.VSA_EState5 = VSA_EState5
        self.VSA_EState6 = VSA_EState6
        self.VSA_EState7 = VSA_EState7
        self.VSA_EState8 = VSA_EState8
        self.VSA_EState9 = VSA_EState9
        self.fr_Al_COO = fr_Al_COO
        self.fr_Al_OH = fr_Al_OH
        self.fr_Al_OH_noTert = fr_Al_OH_noTert
        self.fr_ArN = fr_ArN
        self.fr_Ar_COO = fr_Ar_COO
        self.fr_Ar_N = fr_Ar_N
        self.fr_Ar_NH = fr_Ar_NH
        self.fr_Ar_OH = fr_Ar_OH
        self.fr_COO = fr_COO
        self.fr_COO2 = fr_COO2
        self.fr_C_O = fr_C_O
        self.fr_C_O_noCOO = fr_C_O_noCOO
        self.fr_C_S = fr_C_S
        self.fr_HOCCN = fr_HOCCN
        self.fr_Imine = fr_Imine
        self.fr_NH0 = fr_NH0
        self.fr_NH1 = fr_NH1
        self.fr_NH2 = fr_NH2
        self.fr_N_O = fr_N_O
        self.fr_Ndealkylation1 = fr_Ndealkylation1
        self.fr_Ndealkylation2 = fr_Ndealkylation2
        self.fr_Nhpyrrole = fr_Nhpyrrole
        self.fr_SH = fr_SH
        self.fr_aldehyde = fr_aldehyde
        self.fr_alkyl_carbamate = fr_alkyl_carbamate
        self.fr_alkyl_halide = fr_alkyl_halide
        self.fr_allylic_oxid = fr_allylic_oxid
        self.fr_amide = fr_amide
        self.fr_amidine = fr_amidine
        self.fr_aniline = fr_aniline
        self.fr_aryl_methyl = fr_aryl_methyl
        self.fr_azide = fr_azide
        self.fr_azo = fr_azo
        self.fr_barbitur = fr_barbitur
        self.fr_benzene = fr_benzene
        self.fr_benzodiazepine = fr_benzodiazepine
        self.fr_bicyclic = fr_bicyclic
        self.fr_diazo = fr_diazo
        self.fr_dihydropyridine = fr_dihydropyridine
        self.fr_epoxide = fr_epoxide
        self.fr_ester = fr_ester
        self.fr_ether = fr_ether
        self.fr_furan = fr_furan
        self.fr_guanido = fr_guanido
        self.fr_halogen = fr_halogen
        self.fr_hdrzine = fr_hdrzine
        self.fr_hdrzone = fr_hdrzone
        self.fr_imidazole = fr_imidazole
        self.fr_imide = fr_imide
        self.fr_isocyan = fr_isocyan
        self.fr_isothiocyan = fr_isothiocyan
        self.fr_ketone = fr_ketone
        self.fr_ketone_Topliss = fr_ketone_Topliss
        self.fr_lactam = fr_lactam
        self.fr_lactone = fr_lactone
        self.fr_methoxy = fr_methoxy
        self.fr_morpholine = fr_morpholine
        self.fr_nitrile = fr_nitrile
        self.fr_nitro = fr_nitro
        self.fr_nitro_arom = fr_nitro_arom
        self.fr_nitro_arom_nonortho = fr_nitro_arom_nonortho
        self.fr_nitroso = fr_nitroso
        self.fr_oxazole = fr_oxazole
        self.fr_oxime = fr_oxime
        self.fr_para_hydroxylation = fr_para_hydroxylation
        self.fr_phenol = fr_phenol
        self.fr_phenol_noOrthoHbond = fr_phenol_noOrthoHbond
        self.fr_phos_acid = fr_phos_acid
        self.fr_phos_ester = fr_phos_ester
        self.fr_piperdine = fr_piperdine
        self.fr_piperzine = fr_piperzine
        self.fr_priamide = fr_priamide
        self.fr_prisulfonamd = fr_prisulfonamd
        self.fr_pyridine = fr_pyridine
        self.fr_quatN = fr_quatN
        self.fr_sulfide = fr_sulfide
        self.fr_sulfonamd = fr_sulfonamd
        self.fr_sulfone = fr_sulfone
        self.fr_term_acetylene = fr_term_acetylene
        self.fr_tetrazole = fr_tetrazole
        self.fr_thiazole = fr_thiazole
        self.fr_thiocyan = fr_thiocyan
        self.fr_thiophene = fr_thiophene
        self.fr_unbrch_alkane = fr_unbrch_alkane
        self.fr_urea = fr_urea
        self.qed = qed
        self.tags = tags

    def __repr__(self: TMoleculeModel):
        return "<MoleculeModel(id='%s', inchi_key='%s', inchi='%s', smiles='%s')>" % (
            self.id,
            self.inchi_key,
            self.inchi,
            self.smiles,
        )
