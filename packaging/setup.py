#!/usr/bin/env python3
"""
OncoSieve -- pip-installable package setup.

Install in development mode:
    pip install -e packaging/

This registers the `oncosieve` command-line entry point.
"""

from setuptools import setup, find_packages

setup(
    name="oncosieve",
    version="2.1.0",
    description="Pan-Cancer Somatic Variant Whitelist Curation Tool",
    author="Dr Christopher Trethewey",
    author_email="christopher.trethewey@nhs.net",
    url="https://github.com/Trethewey/oncosieve",
    license="AGPL-3.0-or-later",
    python_requires=">=3.10",
    packages=find_packages(where=".."),
    package_dir={"": ".."},
    install_requires=[
        "pandas>=2.0",
        "pyyaml>=6.0",
        "requests>=2.31",
        "pysam>=0.21",
        "openpyxl>=3.1",
        "xlrd>=2.0",
        "polars>=0.20",
        "pyarrow>=14.0",
        "plotly>=5.0",
        "pyliftover>=0.4",
    ],
    entry_points={
        "console_scripts": [
            "oncosieve=build_whitelist:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
)
