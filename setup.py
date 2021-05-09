# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='iDS_Longterm_Emissions',
    version='0.1.0',
    description='Package for simulating long term emissions from iDS waste bodies',
    long_description=readme,
    author='Timo Heimovaara',
    author_email='t.j.heimovaara@tudelft.nl',
    url='',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

