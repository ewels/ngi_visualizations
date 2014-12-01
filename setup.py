from setuptools import setup, find_packages
import sys, os

version = '0.1'

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(name='ngi_visualizations',
    version=version,
    description="Data visualizations for NGI Sweden",
    long_description='This package is used to generate data visualizations '
                   'of various inputs for the NGI Sweden.',
    keywords='bioinformatics',
    author='Phil Ewels',
    author_email='phil.ewels@scilifelab.se',
    url='https://portal.scilifelab.se/genomics/',
    license='MIT',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requires
)