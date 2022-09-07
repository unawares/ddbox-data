#!/bin/bash

case "$ENV" in
"DEV")
    uvicorn app.main:app --app-dir=./src --host 0.0.0.0 --port 8000 --reload
    ;;
"PRODUCTION")
    uvicorn app.main:app --app-dir=./src --host 0.0.0.0 --port 8000 --proxy-headers --workers 5 # workers = (2*CPU)+1
    ;;
*)
    echo "NO ENV SPECIFIED!"
    exit 1
    ;;
esac
