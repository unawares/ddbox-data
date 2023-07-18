FROM ubuntu:22.04 AS builder

RUN apt-get update \
 && apt-get install -yq --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    wget \
    libboost-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-serialization-dev \
    libboost-python-dev \
    libboost-regex-dev \
    libboost-iostreams-dev \
    libcairo2-dev \
    libeigen3-dev \
    software-properties-common \
    swig \
    python3.10 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://bootstrap.pypa.io/get-pip.py \
 && python3.10 get-pip.py \
 && rm get-pip.py


# Build Open Babel in python3
WORKDIR /
ARG OPENBABEL_VERSION=openbabel-3-1-1
RUN wget --quiet https://github.com/openbabel/openbabel/archive/${OPENBABEL_VERSION}.tar.gz \
 && tar -xzf ${OPENBABEL_VERSION}.tar.gz \
 && mv openbabel-${OPENBABEL_VERSION} openbabel \
 && rm ${OPENBABEL_VERSION}.tar.gz

RUN mkdir /openbabel/build
WORKDIR /openbabel/build

RUN cmake \
  -D PYTHON_BINDINGS=ON \
  -D RUN_SWIG=ON \
  -D PYTHON_EXECUTABLE=/usr/bin/python3.10 \
  -D PYTHON_INCLUDE_DIR=/usr/include/python3.10 \
  -D CMAKE_INSTALL_PREFIX=/usr \
  -D CMAKE_BUILD_TYPE=Release \
  ..

RUN make -j $(nproc)

RUN make install

RUN ln -s /usr/include/openbabel3 /usr/local/include/openbabel3

RUN python3.10 -m pip install openbabel


## Base setup
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONFAULTHANDLER 1
ENV PYTHONUNBUFFERED 1
ENV ROOT /app

ENV PYTHONPATH "${PYTHONPATH}:/usr/lib/python3.10/dist-packages:/app/src"

WORKDIR $ROOT

RUN apt-get update \
  && apt-get install -y --no-install-recommends apt-utils \
  && apt-get install -y --no-install-recommends libc-dev \
  && apt-get install -y --no-install-recommends gcc \
  && apt-get install -y --no-install-recommends gettext \
  && apt-get install -y --no-install-recommends screen \
  && apt-get install -y --no-install-recommends vim \
  && apt-get install -y --no-install-recommends curl \
  && apt-get install -y --no-install-recommends autodock-vina \
  && apt-get clean

COPY requirements.txt ./

RUN python3.10 -m pip install -r requirements.txt

# Add commands
COPY devops/commands $ROOT/commands
RUN chmod +x $ROOT/commands/*
ENV PATH="$ROOT/commands:$PATH"

ADD src $ROOT/src

WORKDIR $ROOT/src

ENTRYPOINT [ "entrypoint.sh" ]
CMD [ "start.sh" ]
