"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='scipion-chem-lanlcatnap',
    version='0.1.0',
    description='Scipion framework plugin cross-referencing peptides against known HIV bnAb epitopes (LANL/CATNAP)',
    long_description=long_description,
    url='https://github.com/Lvera-code/scipion-chem-lanlcatnap',
    author='Enzo Sierra',
    author_email='enzogael57@gmail.com',
    keywords='scipion epitope hiv bnab lanl catnap',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    entry_points={
        'pyworkflow.plugin': 'lanlcatnap = lanlcatnap'
    }
)
