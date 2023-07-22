# parser.py
import csv
import logging
from typing import *

import pandas as pd
from app.db import Session
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

    data = pd.read_csv('../../dataset/data.csv', header=None)

    data = data[STARTINDEX:]

    size = len(data)

    rows: List[List] = []

    BATCH_SIZE = 1000

    def append():
        nonlocal rows

        session = Session()

        tags = {
            'moses': get_or_create(session, TagModel, name='moses')
        }

        for i in range(len(rows)):
            row = rows[i]
            if row[0] not in tags:
                tags[row[0]] = get_or_create(session, TagModel, name=row[0])

            tag_objs = [tags.get(row[0]), tags.get('moses')]

            molecule_obj = MoleculeModel(
                inchi_key=row[1],
                inchi=row[2],
                smiles=row[3],
                **{
                    key: row[4 + index] for index, key in enumerate(MOLECULE_DESCRIPTORS)
                },
                tags=tag_objs,
            )
            session.add(molecule_obj)

        session.commit()
        session.close()

        rows = []

    for index, row in data.iterrows():
        row = row.tolist()
        rows.append(row)
        if len(rows) == BATCH_SIZE:
            append()
            logger.info("Uploaded: %s/%s" % (index + 1, size))

    append()
    logger.info("Uploaded: %s/%s" % (size, size))
