#!/bin/bash

# Usage: bootstrap.sh [virtualenv-name]
#
# Create a virtualenv and install requirements to it
#
set -e

if [[ -z $1 ]]; then
    venv=$(basename $(pwd))-env
else
    venv=$1
fi
mkdir -p src

# Create a virtualenv using a specified version of the virtualenv
# source.  This also provides setuptools and pip.  

# Create virtualenv if necessary
if [ ! -f $venv/bin/activate ]; then
    virtualenv --python /usr/local/bin/python ${venv}
    virtualenv --python /usr/local/bin/python --relocatable ${venv}
else
    echo "found existing virtualenv $venv"
fi

source $venv/bin/activate

# full path; set by activate
venv=$VIRTUAL_ENV

echo $venv
$venv/bin/pip install -r requirements.txt

python setup.py install
