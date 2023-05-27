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


settings = Settings()
