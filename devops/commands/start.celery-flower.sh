#!/bin/bash

celery -A app.celery --broker=$CELERY_BROKER_URL flower -l $CELERY_LOG_LEVEL --port=5555 --basic_auth=$CELERY_FLOWER_USER:$CELERY_FLOWER_PASSWORD
