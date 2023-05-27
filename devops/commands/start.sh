#!/bin/bash

python main.py migrate

case "$ENV" in
"DEV")
    uvicorn main:app --host 0.0.0.0 --port 8000 --reload
    ;;
"PRODUCTION")
    gunicorn main:app --worker-class uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000 --workers $GUNICORN_WORKERS --timeout $GUNICORN_TIMEOUT
    ;;
*)
    echo "NO ENV SPECIFIED!"
    exit 1
    ;;
esac
