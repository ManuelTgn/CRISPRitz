from setuptools import setup, find_packages

import os

# Read README for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [
        line.strip() for line in fh if line.strip() and not line.startswith("#")
    ]

setup(
    name="crispritz",
    version="3.0.0",
    author="Pinellolab",
    author_email="lpinello@mgh.harvard.edu",
    description="CRISPRitz: High-Throughput and Variant-Aware In Silico "
    "Off-Target Sites Identification For CRISPR Genome Editing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pinellolab/CRISPRitz",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": ["pytest", "black", "flake8", "mypy"],
        "docs": ["sphinx", "sphinx-rtd-theme"],
    },
    entry_points={
        "console_scripts": [
            "crispritz=crispritz.__main__:main",
        ],
    },
)
