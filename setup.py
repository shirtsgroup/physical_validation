r"""
`physicalvalidation` is a package aimed at testing results obtained
by molecular dynamics simulations for their physical validity.

shirtsgroup.github.io/physical-validation
"""
from __future__ import print_function

import os
import sys
from setuptools import setup, find_packages

#####################################
VERSION = "0.1.0"
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
    name='physicalvalidation',
    version=__version__,
    description=__doc__,
    author='Michael R. Shirts, Pascal T. Merz',
    author_email='michael.shirts@colorado.edu, pascal.merz@colorado.edu',
    url='https://github.com/shirtsgroup/physical-validation',
    download_url='https://github.com/shirtsgroup/physical-validation/tarball/{}'.format(__version__),
    packages=find_packages(),
    package_dir={'physicalvalidation': 'physicalvalidation'},
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
    license="LGPLv2.1",
    zip_safe=False,
    keywords='physicalvalidation',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: LGPLv2.1 License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Molecular-Simulation',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
)
