# parser.py
import csv

import pandas as pd
from ddbox.molecule.objects import MOLECULE_DESCRIPTORS, Molecule
from rdkit import RDLogger
from tqdm import tqdm

lg = RDLogger.logger()

lg.setLevel(RDLogger.CRITICAL)

data = pd.read_csv('../../../dataset/dataset_v1.csv')

smiles_list = data['SMILES'].to_list()[0:]
split_list = data['SPLIT'].to_list()[0:]

molecules = []

descriptors = []

BATCH_SIZE = 1000


def append():
    global descriptors
    with open('../../../dataset/data.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerows(descriptors)
    descriptors = []


for i in tqdm(range(len(smiles_list))):
    smiles = smiles_list[i]
    split = split_list[i]

    molecule = Molecule.from_smiles(smiles)
    molecules.append(molecule)

    row = [split, molecule.inchi_key, molecule.inchi, molecule.smiles, *[getattr(molecule.desciptor, descriptor) for descriptor in MOLECULE_DESCRIPTORS]]
    descriptors.append(row)

    if len(descriptors) == BATCH_SIZE:
        append()

append()
