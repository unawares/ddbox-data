import logging
import os
from subprocess import PIPE, Popen
from typing import List, Tuple

from app.configs import settings
from services.s3 import S3ServiceBuilder
from utils.dirs import delete_if_exists
from utils.random import get_random_uuid_hex

logger = logging.getLogger(__name__)


def generate_pdbqt_file_from_smiles(smiles: str):
    filepath = '/tmp/%s.pdbqt' % get_random_uuid_hex()
    command = [
        'obabel',
        '-:%s' % smiles,
        '-p', '7.4',
        '--gen3d',
        '-o', 'pdbqt',
        '-O', filepath,
    ]

    p = Popen(command, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    if not ('1 molecule converted' in output.decode() or '1 molecule converted' in err.decode()):
        raise Exception("Invalid smiles")
    return filepath


def generate_pbdqt_file_from_receptor_id(receptor_id: str):
    filename = '%s_target.pdbqt' % receptor_id
    filepath = '/tmp/%s' % filename
    if not (os.path.exists(filepath) and os.path.isfile(filepath)):
        delete_if_exists(filepath)
        s3_service = S3ServiceBuilder().build()
        receptor_content = s3_service.bucket(settings.TARGETS_BUCKET_NAME).download_binary(filename).decode(settings.ENCODING)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(receptor_content)
    return filepath


def vina_docking(receptor_path: str, ligand_path: str, center: Tuple[float, float, float], size: Tuple[float, float, float]):
    command = [
        'vina',
        '--receptor', receptor_path,
        '--ligand', ligand_path,
        '--center_x', str(center[0]),
        '--center_y', str(center[1]),
        '--center_z', str(center[2]),
        '--size_x', str(size[0]),
        '--size_y', str(size[1]),
        '--size_z', str(size[2]),
    ]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    if p.returncode != 0:
        raise Exception(err.decode())
    name, ext = ligand_path.split('.')
    return '%s_out.%s' % (name, ext)
