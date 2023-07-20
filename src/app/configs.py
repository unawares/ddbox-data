import logging

from pydantic import BaseSettings
from sqlalchemy import BigInteger

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(threadName)s] [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)


class Settings(BaseSettings):

    ENV: str
    POSTGRES_HOST: str
    POSTGRES_PORT: str
    POSTGRES_DB: str
    POSTGRES_USER: str
    POSTGRES_PASSWORD: str
    DATABASE_URL: str

    SQLALCHEMY_ID_TYPE = BigInteger

    CELERY_BROKER_URL: str
    CELERY_RESULT_BACKEND: str

    MINIO_ENDPOINT: str
    MINIO_ACCESS_KEY: str
    MINIO_SECRET_KEY: str

    REDIS_CACHE: str

    TARGETS_BUCKET_NAME = 'targets'
    SUBMISSIONS_BUCKET_NAME = 'submissions'
    SUBMISSION_RESULTS_BUCKET_NAME = 'submission-results'
    SUBMISSION_DOCKING_RESULTS_BUCKET_NAME = 'submission-docking-results'
    ENCODING = 'utf-8'


settings = Settings()
