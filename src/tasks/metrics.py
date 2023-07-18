import json
import logging
from typing import *
from typing import List

from app.celery_base import app
from app.configs import settings
from app.db import Session
from celery import shared_task
from ddbox.metrics import metrics_regstry
from ddbox.metrics.utils import get_molecule_from_smiles_if_valid_or_none
from models import MoleculeModel, TagModel
from services.s3 import S3ServiceBuilder

logger = logging.getLogger(__name__)


@app.task
def compute_metrics(submission_id: str, metrics: List[str]):
    s3_service = S3ServiceBuilder().build()

    s3_service \
        .bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME) \
        .upload_binary(
            '%s.json' % submission_id,
            json.dumps({
                'status': 'processing'
            }).encode(settings.ENCODING)
        )

    filepath = '%s.json' % submission_id

    file_content = s3_service.bucket(settings.SUBMISSIONS_BUCKET_NAME).download_binary(filepath).decode(settings.ENCODING)
    smiles_list = json.loads(file_content)

    session = Session()

    tag = session.query(TagModel).filter(TagModel.name == 'test').one()
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

    filepath = '%s.json' % submission_id
    file_content = json.dumps(results).encode(settings.ENCODING)
    s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).upload_binary(filepath, file_content)
