# parser.py
import csv

from ddbox.molecule.descriptors import MOLECULE_DESCRIPTORS
from ddbox.molecule.objects import Molecule
from rdkit import RDLogger
from tools.data.moses import moses_data
from tqdm import tqdm

lg = RDLogger.logger()

lg.setLevel(RDLogger.CRITICAL)

smiles_list = moses_data().pd['SMILES'].to_list()[1792000:]

molecules = []

keys = []
descriptors = []

BATCH_SIZE = 1000


def append():
    global keys, descriptors
    with open('keys.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerows(keys)
    with open('descriptors.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerows(descriptors)
    keys = []
    descriptors = []


for smiles in tqdm(smiles_list):
    molecule = Molecule().from_smiles(smiles)
    molecules.append(molecule)

    keys.append([molecule.inchi_key, molecule.inchi, molecule.smiles])
    descriptors.append(list([getattr(molecule.desciptor, descriptor) for descriptor in MOLECULE_DESCRIPTORS]))

    if len(keys) == BATCH_SIZE:
        append()

append()
