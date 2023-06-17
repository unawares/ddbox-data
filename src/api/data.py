import logging

from app.db import Session
from app.fastapi import fastapi
from models import MoleculeModel, TagModel
from utils.responses import SuccessResponce

logger = logging.getLogger(__name__)


@fastapi.get("/data/moses")
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
