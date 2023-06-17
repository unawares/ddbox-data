#!/bin/bash

celery -A app.celery worker -l $CELERY_LOG_LEVEL -P gevent -c $CELERY_CONCURENCY
