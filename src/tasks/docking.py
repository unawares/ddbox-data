import json
import logging
from typing import List, Tuple

from app.celery_base import app
from app.configs import settings
from ddbox.docking.vina import generate_pbdqt_file_from_receptor_id, generate_pdbqt_file_from_smiles, vina_docking
from services.s3 import S3ServiceBuilder
from utils.dirs import delete_if_exists

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

    docking_results = []
    ligands = {}
    receptors = {}

    ligand_filepaths = {}
    receptor_filepaths = {}

    for i in range(len(receptor_ids)):
        receptor_filepath = generate_pbdqt_file_from_receptor_id(receptor_ids[i])
        receptor_filepaths[receptor_ids[i]] = receptor_filepath

        center = centers[i]
        size = sizes[i]

        for smiles in smiles_list:
            logger.info("Docking %s to receptor %s" % (smiles, receptor_ids[i]))

            if smiles not in ligand_filepaths:
                ligand_filepath = generate_pdbqt_file_from_smiles(smiles)
                ligand_filepaths[smiles] = ligand_filepath
            else:
                ligand_filepath = ligand_filepaths[smiles]

            result_filepath = vina_docking(
                receptor_filepath,
                ligand_filepath,
                center,
                size,
            )

            with open(result_filepath, 'r', encoding='utf-8') as f:
                vina_pdbqt = f.read()

            docking_results.append({
                'smiles': smiles,
                'receptor_id': receptor_ids[i],
                'vina_pdbqt': vina_pdbqt,
            })

    for ligand_smiles in ligand_filepaths.keys():
        with open(ligand_filepaths[smiles], 'r', encoding='utf-8') as f:
            ligands[ligand_smiles] = f.read()

    for receptor_id in receptor_filepaths.keys():
        with open(receptor_filepaths[receptor_id], 'r', encoding='utf-8') as f:
            receptors[receptor_id] = f.read()

    filepath = '%s.json' % submission_id
    file_content = json.dumps({
        'ligands': [
            {
                'smiles': ligand,
                'pdbqt': ligands[ligand]
            } for ligand in ligands
        ],
        'receptors': [
            {
                'receptor_id': receptor,
                'pdbqt': receptors[receptor]
            } for receptor in receptors
        ],
        'docking_results': docking_results
    }).encode(settings.ENCODING)
    s3_service.bucket(settings.SUBMISSION_RESULTS_BUCKET_NAME).upload_binary(filepath, file_content)
