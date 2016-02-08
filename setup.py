"""
Create unix package:    python setup.py sdist
"""

import os
import subprocess
import shutil
from setuptools import setup
from os.path import join

subprocess.call('git log --pretty=format:%h -n 1 > munging/data/sha', shell = True)
subprocess.call('git shortlog --format="XXYYXX%h" | grep -c XXYYXX > munging/data/ver', shell = True)

from munging import __version__
from munging.scripts import script

params = {'author': 'Sheena Scroggins',
          'author_email': 'sheena.scroggins@gmail.com',
          'description': script.__doc__.strip(),
          'name': 'munging',
          'packages': ['munging','munging.scripts','munging.subcommands'],
          'package_dir': {'munging': 'munging'},
          'scripts': ['munge'],
          'version': __version__,
          'package_data': {'munging': [join('data',f) for f in ['sha','ver']]},
          'install_requires': [
              'numpy==1.10.1',
              'pandas==0.17.1',
              'wsgiref==0.1.2',
              'xlrd==0.9.3',
              'xlwt==0.7.5'
          ]
      }

setup(**params)

