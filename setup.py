import os
from setuptools import setup
import glob

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# scriptlist = glob.glob(os.path.join('bin', '*.py'))

setup(
    name="cplutils",
    version="0.0.4",
    author="Eduardo Ramos Fernandez",
    author_email="eduradical951@gmail.com",
    description=("A python package with utilities to do MD-CFD coupled simulations"),
    license="BSD",
    keywords = "multiscale md cfd coupling",
    url = "",
    packages=['cplutils', 'cplutils.readers', 'cplutils.plotting',
              'cplutils.builders', 'cplutils.postproc', 'cplutils.models',
              'cplutils.builders.surfacebuilder'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
    # scripts=scriptlist,
    zip_safe=False,
)
