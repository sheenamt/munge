"""
Test the make_cnv_plottable script
"""

import subprocess
import filecmp
import logging
import os
from munging.subcommands import make_cnv_plottable

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)