#!/bin/bash

# Usage: create_venv.sh <virtualenv> <bioy branch, tag, or commit>

set -e

venv=${1-$(basename $(pwd))-env}
bioy_version=${2-master}

mkdir -p src
if [[ ! -d src/bioy ]]; then
    git clone https://github.com/nhoffman/bioy.git src/bioy
fi

(cd src/bioy && git checkout $bioy_version)

src/bioy/dev/bootstrap.sh \
    --venv $venv \
    --wheelstreet /usr/local/share/python/wheels \
    --requirements src/bioy/requirements.txt

src/bioy/dev/bootstrap.sh \
    --venv $venv \
    --wheelstreet /usr/local/share/python/wheels \
    --requirements requirements.txt

# $venv/bin/pip install -r requirements.txt
