import json
import logging
from typing import List, Tuple

from app.configs import settings
from app.db import Session
from app.fastapi import fastapi
from ddbox.metrics import metrics_regstry
from ddbox.metrics.utils import get_molecule_from_smiles_if_valid_or_none
from models import MoleculeModel, TagModel
from services.s3 import S3ServiceBuilder
from tasks.metrics import compute_metrics
from tasks.docking import dock
from utils.random import get_next_submission_id
from utils.responses import SuccessResponce

logger = logging.getLogger(__name__)


@fastapi.post("/submission/moses")
async def submission_moses(smiles_list: List[str], metrics: List[str] | None = None):
    s3_service = S3ServiceBuilder().build()

    submission_id = get_next_submission_id()
    filepath = '%s.json' % submission_id
    file_content = json.dumps(smiles_list).encode(settings.ENCODING)

    s3_service.bucket(settings.SUBMISSIONS_BUCKET_NAME).upload_binary(filepath, file_content)

    if metrics is not None:
        metrics = [metric for metric in metrics if metric in metrics_regstry.functions]

    if metrics is not None and len(metrics) == 0:
        metrics = None

    if metrics is None:
        metrics = list(metrics_regstry.functions.keys())

    s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).upload_binary(filepath, json.dumps({
        'status': 'pending'
    }).encode(settings.ENCODING))

    compute_metrics.delay(submission_id, metrics)

    return SuccessResponce({
        'submission_id': submission_id
    }).to_response()


@fastapi.get("/submission/result/moses")
async def submission_result_moses(submission_id: str):
    s3_service = S3ServiceBuilder().build()

    filepath = '%s.json' % submission_id

    file_content = s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).download_binary(filepath).decode(settings.ENCODING)
    data = json.loads(file_content)
    return SuccessResponce(data).to_response()


@fastapi.post("/submission/docking")
async def submission_docking(
    smiles_list: List[str],
    receptor_ids: List[str], 
    centers: List[Tuple[float, float, float]],
    sizes: List[Tuple[float, float, float]],
):
    s3_service = S3ServiceBuilder().build()

    submission_id = get_next_submission_id()
    filepath = '%s.json' % submission_id
    file_content = json.dumps(smiles_list).encode(settings.ENCODING)

    s3_service.bucket(settings.SUBMISSIONS_BUCKET_NAME).upload_binary(filepath, file_content)

    s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).upload_binary(filepath, json.dumps({
        'status': 'pending'
    }).encode(settings.ENCODING))

    dock.delay(submission_id, smiles_list, receptor_ids, centers, sizes)

    return SuccessResponce({
        'submission_id': submission_id
    }).to_response()


@fastapi.get("/submission/result/docking")
async def submission_result_docking(submission_id: str):
    s3_service = S3ServiceBuilder().build()

    filepath = '%s.json' % submission_id

    file_content = s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).download_binary(filepath).decode(settings.ENCODING)
    data = json.loads(file_content)
    return SuccessResponce(data).to_response()
