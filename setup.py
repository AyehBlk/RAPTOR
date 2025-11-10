#!/usr/bin/env python3
"""
RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource
Setup configuration for PyPI distribution
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_file(filename):
    """Read file contents."""
    with open(os.path.join(os.path.dirname(__file__), filename), encoding='utf-8') as f:
        return f.read()

# Read version from __init__.py
def get_version():
    """Extract version from raptor/__init__.py"""
    with open('raptor/__init__.py', 'r') as f:
        for line in f:
            if line.startswith('__version__'):
                return line.split('=')[1].strip().strip('"').strip("'")
    return '2.0.0'  # Fallback version

# Core dependencies (minimal installation)
INSTALL_REQUIRES = [
    'numpy>=1.21.0',
    'pandas>=1.3.0',
    'scipy>=1.7.0',
    'matplotlib>=3.4.0',
    'seaborn>=0.11.0',
    'scikit-learn>=1.0.0',
    'pyyaml>=5.4.0',
    'jinja2>=3.0.0',
    'click>=8.0.0',
    'tqdm>=4.62.0',
    'colorama>=0.4.4',
]

# Optional dependencies for development
EXTRAS_REQUIRE = {
    'dev': [
        'pytest>=7.0.0',
        'pytest-cov>=3.0.0',
        'black>=22.0.0',
        'flake8>=4.0.0',
        'mypy>=0.950',
        'sphinx>=4.0.0',
        'sphinx_rtd_theme>=1.0.0',
    ],
    'docs': [
        'sphinx>=4.0.0',
        'sphinx_rtd_theme>=1.0.0',
        'sphinx-click>=3.0.0',
    ],
}

# All optional dependencies
EXTRAS_REQUIRE['all'] = list(set(sum(EXTRAS_REQUIRE.values(), [])))

setup(
    # Basic package information
    name='raptor-rnaseq',
    version=get_version(),
    description='RNA-seq Analysis Pipeline Testing and Optimization Resource',
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    
    # Author information
    author='Ayeh Bolouki',
    author_email='ayehbolouki1988@gmail.com',
    
    # URLs
    url='https://github.com/AyehBlk/RAPTOR',
    project_urls={
        'Bug Reports': 'https://github.com/AyehBlk/RAPTOR/issues',
        'Source': 'https://github.com/AyehBlk/RAPTOR',
        'Documentation': 'https://github.com/AyehBlk/RAPTOR/tree/main/docs',
    },
    
    # License
    license='MIT',
    
    # Package discovery
    packages=find_packages(exclude=['tests', 'tests.*', 'docs', 'examples']),
    
    # Include non-Python files
    include_package_data=True,
    package_data={
        'raptor': [
            'config/*.yaml',
            'templates/*.html',
            'templates/*.md',
        ],
    },
    
    # Dependencies
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    
    # Python version requirement
    python_requires='>=3.8',
    
    # Entry points (command-line scripts)
    entry_points={
        'console_scripts': [
            'raptor=raptor.cli:main',
        ],
    },
    
    # PyPI classifiers
    classifiers=[
        # Development status
        'Development Status :: 4 - Beta',
        
        # Intended audience
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        
        # Topic
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        
        # License
        'License :: OSI Approved :: MIT License',
        
        # Python versions
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        
        # Operating systems
        'Operating System :: OS Independent',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        
        # Other
        'Natural Language :: English',
    ],
    
    # Keywords for PyPI search
    keywords=[
        'rna-seq',
        'differential-expression',
        'bioinformatics',
        'computational-biology',
        'transcriptomics',
        'pipeline',
        'benchmarking',
        'data-analysis',
        'genomics',
    ],
    
    # Zip safe
    zip_safe=False,
)
