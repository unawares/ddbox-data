import logging
import os

from app.configs import settings
from celery import Celery
from celery.schedules import crontab

logger = logging.getLogger(__name__)

app = Celery('app', include=['tasks.metrics'])
app.conf.broker_url = settings.CELERY_BROKER_URL
app.conf.result_backend = settings.CELERY_RESULT_BACKEND

app.conf.beat_schedule = {
    # 'celery_beat_testing': {
    #     'task': 'tasks.metrics.debug_task',
    #     'schedule': crontab(minute='*/1')
    # }
}

app.conf.timezone = 'UTC'
app.conf.update(task_track_started=True)
