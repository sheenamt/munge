"""
Test the subcommand scripts
"""
import os
from os import path
import unittest
import logging
import pprint
import csv
import sys
import json

from munging.subcommands import create_bed

import __init__ as config
log = logging.getLogger(__name__)

