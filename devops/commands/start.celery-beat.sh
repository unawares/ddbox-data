#!/bin/bash

celery -A app.celery beat -l $CELERY_LOG_LEVEL
