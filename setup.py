#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os
import io

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(os.path.join(here, filename), encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

setup(
    name='chimera_p4',
    version='0.2.3',
    url='https://github.com/josan82/chimera_p4',
    author='José Emilio Sánchez',
    author_email='joseemilio.sanchez@uab.cat',
    description=long_description,
    long_description=long_description,
    packages=find_packages(),
    package_data={
        'chimera_p4': ['*.fdef'],
         },
    include_package_data=True,
    platforms='any',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    install_requires=[],
)
