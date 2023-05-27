import logging
from datetime import date, datetime

from app.configs import settings
from app.db import engine
from models import Base

logger = logging.getLogger(__name__)


def migrate():
    logger.info("RUN: migrate")
    logger.info("DATABASE_URL: %s" % settings.DATABASE_URL)
    Base.metadata.create_all(engine)
