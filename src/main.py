# main.py

import logging
import sys

from api import fastapi
from tools import check_metrics, db, upload_into_db

logger = logging.getLogger(__name__)

app = fastapi

if __name__ == '__main__':
    if len(sys.argv) == 2:
        if sys.argv[1] == 'migrate':
            logger.info("Migrating...")
            db.migrate()
        elif sys.argv[1] == 'upload':
            logger.info("Uploading...")
            upload_into_db.upload()
        elif sys.argv[1] == 'check_metrics':
            logger.info("Checking...")
            check_metrics.check()
