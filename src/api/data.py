import logging

from app.configs import settings
from app.db import Session
from app.fastapi import fastapi
from fastapi_cache.decorator import cache
from models import MoleculeModel, ReceptorModel, TagModel
from services.s3 import S3ServiceBuilder
from utils.responses import SuccessResponce

logger = logging.getLogger(__name__)


@fastapi.get("/data/molecules/info/")
@cache(expire=3600 * 24)
async def data_molecules_info(attributes: str | None = None, tags: str | None = None):
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
        for tag in tags:
            query = query.filter(MoleculeModel.tags.any(TagModel.id.in_([tag.id])))

    allowed_attributes = [column.name for column in MoleculeModel.__table__.columns]
    molecules_count = query.count()

    session.close()

    return SuccessResponce({'total': molecules_count, 'allowed_attributes': allowed_attributes}).to_response()


@fastapi.get("/data/molecules/")
@cache(expire=3600 * 24)
async def data_molecules(attributes: str | None = None, tags: str | None = None, offset: int = 0, limit: int = 100):
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
        for tag in tags:
            query = query.filter(MoleculeModel.tags.any(TagModel.id.in_([tag.id])))

    molecules = query \
        .order_by(MoleculeModel.id) \
        .offset(offset) \
        .limit(limit) \
        .all()

    records = [[getattr(molecule, getattr(MoleculeModel, attribute).name) for attribute in attributes] for molecule in molecules]

    session.close()

    return SuccessResponce({'attributes': attributes, 'records': records}).to_response()


@fastapi.get("/data/docking/receptors/info/")
@cache(expire=3600 * 24)
async def data_docking_receptors_info():
    session = Session()

    query = session.query(ReceptorModel)

    molecules_count = query.count()

    session.close()

    return SuccessResponce({'total': molecules_count}).to_response()


@fastapi.get("/data/docking/receptors/")
@cache(expire=3600 * 24)
async def data_docking_receptors(offset: int = 0, limit: int = 100):
    session = Session()

    query = session.query(ReceptorModel)

    receptors = query \
        .order_by(ReceptorModel.id) \
        .offset(offset) \
        .limit(limit) \
        .all()

    session.close()

    return SuccessResponce({'records': [
        {
            'receptor_id': receptor.receptor_id,
            'vina': {
                'pdbqt_url': receptor.vina_receptor_pdbqt_url,
                'config_url': receptor.vina_receptor_config_url,
            }
        } for receptor in receptors
    ]}).to_response()


@fastapi.get("/data/docking/receptors/{receptor_id}/")
@cache(expire=3600 * 24)
async def data_docking_receptor_by_id(receptor_id: str):
    session = Session()

    receptor = session.query(ReceptorModel).filter_by(receptor_id=receptor_id).one_or_none()

    session.close()

    return SuccessResponce({
        'receptor_id': receptor.receptor_id,
        'vina': {
            'pdbqt_url': receptor.vina_receptor_pdbqt_url,
            'config_url': receptor.vina_receptor_config_url,
        }
    }).to_response()
