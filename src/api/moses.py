import logging

from app.db import Session
from app.fastapi import fastapi
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from models import MoleculeModel
from utils.responses import SuccessResponce

logger = logging.getLogger(__name__)


@fastapi.get("/moses")
async def moses_smiles(attributes: str | None, offset: int = 0, limit: int = 100):
    session = Session()

    if attributes is not None:
        attributes = [attribute.strip() for attribute in attributes.split(',') if hasattr(MoleculeModel, attribute.strip())]

    if len(attributes) == 0:
        attributes = None

    if attributes is None:
        attributes = [column.name for column in MoleculeModel.__table__.columns]

    query = session.query(MoleculeModel).with_entities(*(getattr(MoleculeModel, attribute) for attribute in attributes))

    molecules = query \
        .order_by(MoleculeModel.id) \
        .offset(offset) \
        .limit(limit) \
        .all()

    records = [[getattr(molecule, getattr(MoleculeModel, attribute).name) for attribute in attributes] for molecule in molecules]

    session.close()

    return SuccessResponce({'attributes': attributes, 'records': records}).to_response()
