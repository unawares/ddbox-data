import logging
from typing import List

from app.db import Session
from app.fastapi import fastapi
from ddbox.metrics import metrics_regstry
from ddbox.metrics.utils import get_molecule_from_smiles_if_valid_or_none
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from models import MoleculeModel, TagModel
from utils.responses import SuccessResponce

logger = logging.getLogger(__name__)


@fastapi.get("/moses")
async def moses_smiles(attributes: str | None = None, tags: str | None = None, offset: int = 0, limit: int = 100):
    session = Session()

    if attributes is not None:
        attributes = [attribute.strip() for attribute in attributes.split(',') if hasattr(MoleculeModel, attribute.strip())]

    if tags is not None:
        tags = session.query(TagModel).filter(TagModel.name.in_([tag.strip() for tag in tags.split(',')])).all()

    if attributes is not None and len(attributes) == 0:
        attributes = None

    if attributes is None:
        attributes = [column.name for column in MoleculeModel.__table__.columns]

    query = session.query(MoleculeModel).with_entities(*(getattr(MoleculeModel, attribute) for attribute in attributes))

    if tags is not None and len(tags) == 0:
        tags = None

    if tags is not None:
        query = query.filter(MoleculeModel.tags.any(TagModel.id.in_([tag.id for tag in tags])))

    molecules = query \
        .order_by(MoleculeModel.id) \
        .offset(offset) \
        .limit(limit) \
        .all()

    records = [[getattr(molecule, getattr(MoleculeModel, attribute).name) for attribute in attributes] for molecule in molecules]

    session.close()

    return SuccessResponce({'attributes': attributes, 'records': records}).to_response()


@fastapi.post("/moses/eval")
async def moses_smiles(smiles_list: List[str], metrics: List[str] | None = None):
    session = Session()

    if metrics is not None:
        metrics = [metric for metric in metrics if metric in metrics_regstry.functions]

    tag = session.query(TagModel).filter(TagModel.name == 'test').one()

    if metrics is not None and len(metrics) == 0:
        metrics = None

    if metrics is None:
        metrics = metrics_regstry.functions.keys()

    query = session.query(MoleculeModel).with_entities(MoleculeModel.smiles)
    query = query.filter(MoleculeModel.tags.any(TagModel.id == tag.id))

    test_molecules = query \
        .order_by(MoleculeModel.id) \
        .offset(0) \
        .limit(len(smiles_list)) \
        .all()

    session.close()

    test_molecules = [get_molecule_from_smiles_if_valid_or_none(molecule.smiles) for molecule in test_molecules]

    molecules = [get_molecule_from_smiles_if_valid_or_none(smiles) for smiles in smiles_list]
    molecules = [molecule for molecule in molecules if molecule is not None]

    results = {}

    logger.info("smiles_list size: %s" % len(smiles_list))
    logger.info("test_molecules size: %s" % len(test_molecules))

    if 'fraction_valid' in metrics:
        logger.info("Computing: fraction_valid")
        results['fraction_valid'] = metrics_regstry.eval('fraction_valid', smiles_list)

    if 'fraction_unique' in metrics:
        logger.info("Computing: fraction_unique")
        results['fraction_unique'] = metrics_regstry.eval('fraction_unique', smiles_list)

    if 'fraction_passes' in metrics:
        logger.info("Computing: fraction_passes")
        results['fraction_passes'] = metrics_regstry.eval('fraction_passes', molecules)

    if 'IntDiv' in metrics:
        logger.info("Computing: IntDiv")
        results['IntDiv'] = metrics_regstry.eval('IntDiv', molecules)

    if 'IntDiv2' in metrics:
        logger.info("Computing: IntDiv2")
        results['IntDiv2'] = metrics_regstry.eval('IntDiv2', molecules)

    if 'FCD' in metrics:
        logger.info("Computing: FCD")
        results['FCD'] = metrics_regstry.eval('FCD', molecules, test_molecules)

    if 'SNN' in metrics:
        logger.info("Computing: SNN")
        results['SNN'] = metrics_regstry.eval('SNN', molecules, test_molecules)

    if 'Frag' in metrics:
        logger.info("Computing: Frag")
        results['Frag'] = metrics_regstry.eval('Frag', molecules, test_molecules)

    if 'Scaf' in metrics:
        logger.info("Computing: Scaf")
        results['Scaf'] = metrics_regstry.eval('Scaf', molecules, test_molecules)

    if 'DistributionDifferenceLogP' in metrics:
        logger.info("Computing: DistributionDifferenceLogP")
        results['DistributionDifferenceLogP'] = metrics_regstry.eval('DistributionDifferenceLogP', molecules, test_molecules)

    if 'DistributionDifferenceSA' in metrics:
        logger.info("Computing: DistributionDifferenceSA")
        results['DistributionDifferenceSA'] = metrics_regstry.eval('DistributionDifferenceSA', molecules, test_molecules)

    if 'DistributionDifferenceQED' in metrics:
        logger.info("Computing: DistributionDifferenceQED")
        results['DistributionDifferenceQED'] = metrics_regstry.eval('DistributionDifferenceQED', molecules, test_molecules)

    if 'DistributionDifferenceWeight' in metrics:
        logger.info("Computing: DistributionDifferenceWeight")
        results['DistributionDifferenceWeight'] = metrics_regstry.eval('DistributionDifferenceWeight', molecules, test_molecules)

    return SuccessResponce(results).to_response()
