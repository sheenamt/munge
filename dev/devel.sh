#!/bin/bash

# Usage: devel.sh <virtualenv>

set -e

venv=${1-$(basename $(pwd))-env}

# Install dependencies to it
dev/bootstrap.sh $venv

