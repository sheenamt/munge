"""
Create unix package:    python setup.py sdist
"""

import os
import subprocess
import shutil
from setuptools import setup
from os.path import join
import glob

subprocess.call('git log --pretty=format:%h -n 1 > munging/data/sha', shell = True)
subprocess.call('git shortlog --format="XXYYXX%h" | grep -c XXYYXX > munging/data/ver', shell = True)

from munging import __version__
from munging.scripts import script
package_data=glob.glob('munging/data/*')
params = {'author': 'Sheena Todhunter',
          'author_email': 'sheena.todhunter@gmail.com',
          'description': script.__doc__.strip(),
          'name': 'munging',
          'packages': ['munging','munging.scripts','munging.subcommands'],
          'package_dir': {'munging': 'munging'},
          'scripts': ['munge'],
          'version': __version__,
          'setup_requires': 'numpy>=1.10.1',
          'install_requires': [
              'pandas>=0.17.1',              
              'numpy>=1.10.1',
              'natsort>=5.0.2',
              'suds==0.4',
              'XlsxWriter==0.9.6',
              'xlwt==0.7.5',
              'xlrd==0.9.3',
          ]
      }

setup(**params)

