from datetime import date, datetime
from typing import *

from app.configs import settings
from models.base import Base
from services.s3 import S3ServiceBuilder
from sqlalchemy import URL, BigInteger, Column, Date, DateTime, Float, ForeignKey, Integer, String, Table
from sqlalchemy.orm import relationship

TReceptorModel = TypeVar('TReceptorModel', bound='ReceptorModel')


class ReceptorModel(Base):

    __tablename__ = 'receptors'

    id = Column(settings.SQLALCHEMY_ID_TYPE, primary_key=True, index=True)
    receptor_id = Column(String, unique=True)
    vina_receptor_pdbqt_path = Column(String)
    vina_receptor_config_path = Column(String)

    @property
    def vina_receptor_pdbqt_url(self):
        s3_service = S3ServiceBuilder().build()
        url = s3_service.bucket(settings.TARGETS_BUCKET_NAME).get_presigned_url(self.vina_receptor_pdbqt_path)
        return url

    @property
    def vina_receptor_config_url(self):
        s3_service = S3ServiceBuilder().build()
        url = s3_service.bucket(settings.TARGETS_BUCKET_NAME).get_presigned_url(self.vina_receptor_config_path)
        return url

    def __init__(
            self: TReceptorModel,
            receptor_id: str,
            vina_receptor_pdbqt_path: str,
            vina_receptor_config_path: str,
    ):
        self.receptor_id = receptor_id
        self.vina_receptor_pdbqt_path = vina_receptor_pdbqt_path
        self.vina_receptor_config_path = vina_receptor_config_path

    def __repr__(self: TReceptorModel):
        return "<ReceptorModel(id='%s', receptor_id='%s')>" % (
            self.id,
            self.receptor_id,
        )
