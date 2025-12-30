#!/usr/bin/env python3
"""
RAPTOR - RNA-seq Analysis Pipeline Testing and Optimization Resource

Setup configuration for PyPI distribution.

Version: 2.1.2 (Hotfix - Python 3.8-3.11 compatibility)

Author: Ayeh Bolouki
License: MIT
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

# Version
VERSION = '2.1.2'

# Core dependencies
INSTALL_REQUIRES = [
    'numpy>=1.21.0',
    'pandas>=1.3.0',
    'scipy>=1.7.0',
    'matplotlib>=3.4.0',
    'seaborn>=0.11.0',
    'scikit-learn>=1.0.0',
    'joblib>=1.1.0',
    'pyyaml>=5.4.0',
    'jinja2>=3.0.0',
    'click>=8.0.0',
    'tqdm>=4.62.0',
    'colorama>=0.4.4',
    'statsmodels>=0.13.0',
]

# Optional dependencies
EXTRAS_REQUIRE = {
    'ml': [
        'scikit-learn>=1.0.0',
        'joblib>=1.1.0',
    ],
    'dashboard': [
        'streamlit>=1.28.0',
        'plotly>=5.14.0',
        'psutil>=5.8.0',
    ],
    'advanced': [
        'bayesian-optimization>=1.2.0',
        'markdown>=3.3.0',
        'weasyprint>=52.5',
        'colorlog>=6.6.0',
    ],
    'dev': [
        'pytest>=7.0.0',
        'pytest-cov>=3.0.0',
        'black>=22.0.0',
        'flake8>=4.0.0',
        'mypy>=0.950',
    ],
}

# All optional dependencies
EXTRAS_REQUIRE['all'] = list(set(
    dep for deps in EXTRAS_REQUIRE.values() for dep in deps
))

setup(
    name='raptor-rnaseq',
    version=VERSION,
    author='Ayeh Bolouki',
    author_email='ayeh.bolouki@unamur.be',
    description='RNA-seq Analysis Pipeline Testing and Optimization Resource with ML-powered recommendations, adaptive threshold optimization, and Python 3.8+ compatibility',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/AyehBlk/RAPTOR',
    project_urls={
        'Bug Tracker': 'https://github.com/AyehBlk/RAPTOR/issues',
        'Documentation': 'https://github.com/AyehBlk/RAPTOR/tree/main/docs',
        'Source': 'https://github.com/AyehBlk/RAPTOR',
        'Changelog': 'https://github.com/AyehBlk/RAPTOR/blob/main/docs/CHANGELOG.md',
    },
    packages=find_packages(exclude=['tests', 'tests.*', 'docs', 'examples']),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    entry_points={
        'console_scripts': [
            'raptor=raptor.cli:main',
        ],
    },
    include_package_data=True,
    package_data={
        'raptor': [
            'config/*.yaml',
            'templates/*.html',
            'templates/*.md',
        ],
    },
    keywords=[
        'rna-seq',
        'differential-expression',
        'bioinformatics',
        'transcriptomics',
        'pipeline',
        'benchmarking',
        'machine-learning',
        'threshold-optimization',
    ],
    zip_safe=False,
)
