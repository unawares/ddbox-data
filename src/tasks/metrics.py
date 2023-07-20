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
from ddbox.metrics.benchmarks.moses import compute_metrics as compute_moses_metrics

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

    logger.info("smiles_list size: %s" % len(smiles_list))
    logger.info("test_molecules size: %s" % len(test_molecules))

    results = compute_moses_metrics(smiles_list, [molecule.smiles for molecule in test_molecules], metrics)
    
    filepath = '%s.json' % submission_id
    file_content = json.dumps(results).encode(settings.ENCODING)
    s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).upload_binary(filepath, file_content)
