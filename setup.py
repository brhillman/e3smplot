#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="e3smplot",
    version="0.0.1",
    author="Benjamin R. Hillman",
    author_email="bhillma@sandia.gov",
    description="Tools for plotting earth system model output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/brhillman/e3smplot",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: None",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
