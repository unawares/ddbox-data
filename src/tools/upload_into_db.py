# parser.py
import csv
import logging
from typing import *

from app.db import Session
from data.moses import moses_data
from ddbox.molecule.descriptors import MOLECULE_DESCRIPTORS
from ddbox.molecule.objects import Molecule
from models import MoleculeModel, TagModel
from rdkit import RDLogger
from tqdm import tqdm
from utils.orm import get_or_create

logger = logging.getLogger(__name__)


def upload():
    STARTINDEX = 0

    lg = RDLogger.logger()

    lg.setLevel(RDLogger.CRITICAL)

    md = moses_data()

    smiles_list = md.pd['SMILES'].to_list()[STARTINDEX:]
    split_list = md.pd['SPLIT'].to_list()[STARTINDEX:]

    molecules: List[Molecule] = []
    splits: List[str] = []

    BATCH_SIZE = 1000

    def append():
        nonlocal molecules, splits

        session = Session()

        tags = {}

        for i in range(len(molecules)):
            if splits[i] not in tags:
                tags[splits[i]] = get_or_create(session, TagModel, name=splits[i])

            tag_objs = [tags.get(splits[i])]

            mol = molecules[i]

            molecule_obj = MoleculeModel(
                inchi_key=mol.inchi_key,
                inchi=mol.inchi,
                smiles=mol.smiles,
                **{
                    key: getattr(mol.desciptor, key) for key in MOLECULE_DESCRIPTORS
                },
                tags=tag_objs,
            )
            session.add(molecule_obj)

        session.commit()
        session.close()

        molecules = []
        splits = []

    for index, smiles in enumerate(smiles_list):
        molecule = Molecule().from_smiles(smiles)
        molecules.append(molecule)
        splits.append(split_list[index])

        if len(molecules) == BATCH_SIZE:
            append()
            logger.info("Uploaded: %s/%s" % (index + 1, len(smiles_list)))

    append()
    logger.info("Uploaded: %s/%s" % (len(smiles_list), len(smiles_list)))
