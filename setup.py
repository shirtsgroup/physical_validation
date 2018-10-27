r"""
`physical_validation` is a package aimed at testing results obtained
by molecular simulations for their physical validity.

shirtsgroup.github.io/physical-validation
"""
from __future__ import print_function

from setuptools import setup, find_packages
from os import path

#####################################
VERSION = "1.0.0rc7"
__version__ = VERSION

#####################################

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='physical_validation',
    version=__version__,
    description='Physical validation of molecular simulation results',
    long_description='\n' + long_description,
    url='https://physical-validation.readthedocs.io',
    author='Pascal T. Merz and Michael R. Shirts',
    author_email='pascal.merz@colorado.edu, michael.shirts@colorado.edu',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Environment :: Console',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
    keywords='physical-validation, molecular-simulation, molecular-dynamics, molecular-mechanics',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pymbar'
        ],
    package_dir={'physical_validation': 'physical_validation'},
    project_urls={
        'Bug Reports': 'https://github.com/shirtsgroup/physical_validation/issues',
        'Documentation': 'https://physical-validation.readthedocs.io',
        'Source': 'https://github.com/shirtsgroup/physical_validation',
    },
    license="LGPLv2.1",
)
