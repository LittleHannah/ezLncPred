# -*- coding: utf-8 -*-


"""setup.py: setuptools control."""

import re
from setuptools import setup

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('longdist/longdist.py').read(),
    re.M
).group(1)

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")

setup(
    name="longdist",
    packages=["longdist"],
    entry_points={
        "console_scripts": ['longdist = longdist.longdist:main']
    },
    install_requires=[
        'biopython', 'scipy', 'numpy', 'scikit-learn', 'matplotlib', 'configparser', 'sklearn'
    ],
    version=version,
    description="longdist: Method implementation for long ncRNAs and PCT distinction. This application can create and use models base on the method by Schneider et al (2017).",
    long_description=long_descr,
    author="Hugo Wruck Schneider",
    author_email="hugowschneider@gmail.com",
    url="https://github.com/hugowschneider/longdist.py",
)
