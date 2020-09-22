"""
Setup file for the *Twind* package.

"""


# %% IMPORTS
# Built-in imports
from codecs import open
import re

# Package imports
from setuptools import find_packages, setup


# %% SETUP DEFINITION
# Get the long description from the README file
with open('README.md', 'r') as f:
    long_description = f.read()

# Get the requirements list
with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()

# Read the __version__.py file
with open('twind/__version__.py', 'r') as f:
    vf = f.read()

# Obtain version from read-in __version__.py file
version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", vf, re.M).group(1)

# Setup function declaration
setup(name="twind",
      version=version,
      author="Chang-Goo Kim",
      author_email="changgookim@gmail.com",
      description=("Twind: TIGRESS Wind Launching Model and Sampler"),
      long_description=long_description,
      url="https://twind.readthedocs.io",
      project_urls={
          'Documentation': "https://twind.readthedocs.io",
          'Source Code': "https://github.com/changgoo/Twind",
          },
      license='MIT',
      platforms=['Windows', 'Mac OS-X', 'Linux', 'Unix'],
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          ],
      keywords=("galactic wind"),
      python_requires='>=3.6, <4',
      packages=find_packages(),
      package_dir={'twind': "twind"},
      include_package_data=True,
      install_requires=requirements,
      zip_safe=False,
      )
