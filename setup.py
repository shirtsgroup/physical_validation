"""
physical_validation
A package aimed at testing results obtained by molecular simulations for their physical validity.
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except IOError:
    long_description = "\n".join(short_description[2:]),


setup(
    name='physical_validation',
    author='Pascal T. Merz and Michael R. Shirts',
    author_email='pascal.merz@me.com, michael.shirts@colorado.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    url='https://physical-validation.readthedocs.io',
    install_requires=[
        'numpy',
        'scipy',
        'pymbar',
        'matplotlib'
        ],
    project_urls={
        'Bug Reports': 'https://github.com/shirtsgroup/physical_validation/issues',
        'Documentation': 'https://physical-validation.readthedocs.io',
        'Source': 'https://github.com/shirtsgroup/physical_validation',
    },

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Environment :: Console',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],

    keywords='physical-validation, molecular-simulation, molecular-dynamics, molecular-mechanics',
)
