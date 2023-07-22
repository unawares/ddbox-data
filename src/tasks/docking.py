import json
import logging
import os
from typing import List, Tuple

from app.celery_base import app
from app.configs import settings
from app.db import Session
from ddbox.docking.utils.vina import get_pdbqt_file_from_smiles
from ddbox.docking.vina import get_pdbqt_file_from_smiles, vina_docking_files
from models import ReceptorModel
from services.s3 import S3ServiceBuilder
from utils.dirs import ensure_path, exists_file

logger = logging.getLogger(__name__)


@app.task
def dock(submission_id: str, smiles_list: List[str], receptor_ids: List[str], centers: List[Tuple], sizes: List[Tuple]):
    s3_service = S3ServiceBuilder().build()

    s3_service \
        .bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME) \
        .upload_binary(
            '%s.json' % submission_id,
            json.dumps({
                'status': 'processing'
            }).encode(settings.ENCODING)
        )

    vina_docking_results = []
    ligand_pdbqt_filepaths = {}

    session = Session()

    for i in range(len(receptor_ids)):

        receptor = session.query(ReceptorModel).filter_by(receptor_id=receptor_ids[i]).one_or_none()

        if receptor is not None:
            receptor_pdbqt_filepath = os.path.join(settings.TEMP_DIR, 'receptors/%s_target.pdbqt' % receptor_ids[i])
            ensure_path(os.path.dirname(receptor_pdbqt_filepath))

            if not exists_file(receptor_pdbqt_filepath):
                file_data = s3_service.bucket(settings.TARGETS_BUCKET_NAME) \
                    .download_binary(receptor.vina_receptor_pdbqt_path)
                with open(receptor_pdbqt_filepath, 'wb') as f:
                    f.write(file_data)

            receptor_config_filepath = os.path.join(settings.TEMP_DIR, 'receptors/%s_conf.txt' % receptor_ids[i])
            ensure_path(os.path.dirname(receptor_config_filepath))

            if not exists_file(receptor_config_filepath):
                file_data = s3_service.bucket(settings.TARGETS_BUCKET_NAME) \
                    .download_binary(receptor.vina_receptor_config_path)
                with open(receptor_config_filepath, 'wb') as f:
                    f.write(file_data)

            center = centers[i]
            size = sizes[i]

            for smiles in smiles_list:
                logger.info("Docking %s to receptor %s" % (smiles, receptor_ids[i]))

                if smiles not in ligand_pdbqt_filepaths:
                    ligand_pdbqt_filepath = get_pdbqt_file_from_smiles(smiles, temp_dir=os.path.join(settings.TEMP_DIR, 'vina/'))
                    ligand_pdbqt_filepaths[smiles] = ligand_pdbqt_filepath
                else:
                    ligand_pdbqt_filepath = ligand_pdbqt_filepaths[smiles]

                result_filepath = vina_docking_files(
                    receptor_pdbqt_filepath,
                    ligand_pdbqt_filepath,
                    receptor_config_filepath,
                    center,
                    size,
                )

                with open(result_filepath, 'r', encoding='utf-8') as f:
                    vina_pdbqt = f.read()

                vina_docking_results.append({
                    'smiles': smiles,
                    'receptor_id': receptor_ids[i],
                    'vina_pdbqt': vina_pdbqt,
                })

    filepath = '%s.json' % submission_id
    file_content = json.dumps({
        'results': {
            'vina': {
                'docking': vina_docking_results,
            },
        }
    }).encode(settings.ENCODING)
    s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).upload_binary(filepath, file_content)
