from app.configs import settings
from fastapi import FastAPI

from fastapi_cache import FastAPICache
from fastapi_cache.backends.redis import RedisBackend

from redis import asyncio as aioredis

fastapi = FastAPI()


@fastapi.on_event("startup")
async def startup():
    redis = aioredis.from_url(settings.REDIS_CACHE)
    FastAPICache.init(RedisBackend(redis), prefix='fastapi-cache')
