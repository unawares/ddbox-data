#!/bin/bash

case "$ENV" in
"DEV")
    uvicorn src.app.main:app --host 0.0.0.0 --port 8000 --reload
    ;;
"PRODUCTION")
    uvicorn src.app.main:app --host 0.0.0.0 --port 8000 --proxy-headers --workers 5 # workers = (2*CPU)+1
    ;;
*)
    echo "NO ENV SPECIFIED!"
    exit 1
    ;;
esac
