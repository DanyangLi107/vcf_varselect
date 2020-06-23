
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from codecs import open
import os

loc = os.path.abspath(os.path.dirname(__file__))
bindir = os.path.join(loc, "bin/")

with open(os.path.join(loc, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='vcf_varselect',
    version='1.0.0',
    description='variants selection from vcf file',
    author = 'Danyang Li',
    author_email = 'danyang.li@ki.se',
    long_description = long_description,
    url = '',
    license = 'MIT License',
    packages = ['vcf_varselect'],
    keywords = [
        'Variants'
        'VCF',
        'Whole-exome sequencing'
    ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Operating System :: Windows",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)