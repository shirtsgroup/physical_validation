"""checkensemble: A tool for validating molecular simulations.
   See "Simple Quantitative Tests to Validate Sampling from Thermodynamic Ensembles"
   http://pubs.acs.org/doi/abs/10.1021/ct300688p
"""
from __future__ import print_function

import os
import sys
from setuptools import setup, find_packages

#####################################
VERSION = "0.9.0"
ISRELEASED = False
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################

with open('intermol/version.py', 'w') as version_file:
    version_file.write('version="{0}"\n'.format(__version__))

with open('__conda_version__.txt', 'w') as conda_version:
    conda_version.write(__version__)

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

with open('requirements.txt') as reqs_file:
    reqs = [line.strip() for line in reqs_file]

setup(
    name='checkensemble',
    version=__version__,
    description=__doc__,
    author='Michael Shirts, Pascal Merz',
    author_email='michael.shirts@colorado.edu'
    url='https://github.com/shirtsgroup/checkensemble',
    download_url='https://github.com/shirtsgroup/checkensemble/tarball/{}'.format(__version__),
    packages=find_packages(),
    package_dir={'checkensemble': 'checkensemble'},
    package_data={'tests': ['*.py',
                            '*.md',
                            'desmond/*.cfg',
                            'desmond/*/*.cms',
                            'gromacs/*.mdp',
                            'gromacs/*/*/*.gro',
                            'gromacs/*/*/*.top',
                            'gromacs/*/*/*.itp',
                            'lammps/*/*.lmp',
                            'lammps/*/*.input',
                            ]},
    include_package_data=True,
    data_files=[('my_data', ['data/data_file'])],
    install_requires=reqs,
    license="GPL",
    zip_safe=False,
    keywords='checkensemble',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
)

