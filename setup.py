from setuptools import setup
from pathlib import Path

def get_install_requires():
    """Returns requirements.txt parsed to a list"""
    fname = Path(__file__).parent / 'requires.txt'
    targets = []
    if fname.exists():
        with open(fname, 'r') as f:
            targets = f.read().splitlines()
    return targets

setup(name='HiPGx',
      version='0.1.1',
      description='Star Allele Caller for HiFi sequencing',
      author='John Harting',
      author_email='jharting@pacb.com',
      url='http://bitbucket.pacificbiosciences.com:7990/users/jharting/repos/hifi_pgx',
      install_requires=get_install_requires(),
     )
