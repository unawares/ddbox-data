ARG PYTHON_VERSION

FROM python:$PYTHON_VERSION AS base

## Base setup
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONFAULTHANDLER 1
ENV PYTHONBUFFERED 1
ENV ROOT /app

RUN apt-get update \
  && apt-get install -y --no-install-recommends apt-utils \
  && apt-get install -y --no-install-recommends libc-dev \
  && apt-get install -y --no-install-recommends gcc \
  && apt-get install -y --no-install-recommends gettext \
  && apt-get install -y --no-install-recommends screen \
  && apt-get clean

RUN pip install --upgrade pip

## Install pip dependencies
COPY requirements.txt $ROOT/requirements.txt
RUN pip install --no-cache-dir --compile -Ur $ROOT/requirements.txt

# Add commands
COPY devops/commands $ROOT/commands
RUN chmod +x $ROOT/commands/*
ENV PATH="$ROOT/commands:$PATH"

ADD src $ROOT/src

WORKDIR $ROOT

ENTRYPOINT [ "entrypoint.sh" ]
