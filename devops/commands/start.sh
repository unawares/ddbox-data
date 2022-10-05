#!/bin/bash

case "$ENV" in
"DEV")
    uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
    ;;
"PRODUCTION")
    gunicorn app.main:app --worker-class uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000 --workers 5  # workers = (2*CPU)+1
    ;;
*)
    echo "NO ENV SPECIFIED!"
    exit 1
    ;;
esac
