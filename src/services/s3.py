import io
import logging
from dataclasses import dataclass
from datetime import timedelta

import urllib3
from app.configs import settings
from minio import Minio

logger = logging.getLogger(__name__)


@dataclass
class S3Service:
    client: Minio
    bucket_name: str = None

    def bucket(self, bucket_name: str):
        self.bucket_name = bucket_name

        if not self.client.bucket_exists(self.bucket_name):
            self.client.make_bucket(self.bucket_name)

        return self

    def upload_binary(self, filepath: str, binary: bytes):
        stream = io.BytesIO(binary)
        self.client.put_object(self.bucket_name,
                               filepath,
                               stream,
                               length=len(binary))

    def download_binary(self, filepath: str) -> bytes:
        response = None
        try:
            response = self.client.get_object(self.bucket_name, filepath)
            return response.read()
        finally:
            if response is not None:
                response.close()
                response.release_conn()

    def get_presigned_url(self, filepath: str, expiry_seconds: int = 86400) -> str:
        expires = timedelta(seconds=expiry_seconds)
        return self.client.get_presigned_url(method='GET',
                                             bucket_name=self.bucket_name,
                                             object_name=filepath,
                                             expires=expires)


@dataclass
class S3ServiceBuilder:
    host: str = settings.MINIO_ENDPOINT
    access_key: str = settings.MINIO_ACCESS_KEY
    secret_key: str = settings.MINIO_SECRET_KEY

    def build(self):
        # TODO: create try catch for S3Error
        client = Minio(
            self.host,
            access_key=self.access_key,
            secret_key=self.secret_key,
            secure=False,
        )

        return S3Service(client)
