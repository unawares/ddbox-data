from datetime import date, datetime
from typing import *

from app.configs import settings
from models.base import Base
from sqlalchemy import URL, BigInteger, Column, Date, DateTime, Float, ForeignKey, Integer, String, Table
from sqlalchemy.orm import relationship

# TTSmilesListSubmission = TypeVar('TTagModel', bound='TagModel')


# class TagModel(Base):
#     __tablename__ = 'tags'

#     id = Column(settings.SQLALCHEMY_ID_TYPE, primary_key=True, index=True)
#     name = Column(String, unique=True)

#     def __init__(
#             self: TTagModel,
#             name: str,
#     ):
#         self.name = name

#     def __repr__(self: TTagModel):
#         return "<TagModel(id='%s', name='%s')>" % (
#             self.id,
#             self.name,
#         )


# TODO: SmilesListSubmission -> { id, submission_id, json_file }
# TODO: SmilesListSubmissionEvaluation -> { id, metric, value, status: PENDING, PROCESSING, FINISHED }
