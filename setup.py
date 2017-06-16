r"""
`physical_validation` is a package aimed at testing results obtained
by molecular dynamics simulations for their physical validity.

shirtsgroup.github.io/physical-validation
"""
from __future__ import print_function

from setuptools import setup, find_packages

#####################################
VERSION = "0.1.0"
ISRELEASED = False
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + 'a1'
#####################################

with open('requirements.txt') as reqs_file:
    reqs = [line.strip() for line in reqs_file]

setup(
    name='physical_validation',
    version=__version__,
    description=__doc__,
    author='Michael R. Shirts, Pascal T. Merz',
    author_email='michael.shirts@colorado.edu, pascal.merz@colorado.edu',
    url='https://github.com/shirtsgroup/physical-validation',
    packages=find_packages(),
    package_dir={'physical_validation': 'physical_validation'},
    install_requires=reqs,
    license="LGPLv2.1",
    zip_safe=False,
    keywords='physical_validation',
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
