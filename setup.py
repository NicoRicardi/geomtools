#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:28:40 2019

@author: nico
"""
from distutils.core import setup
setup(
    name='geomtools',
    version='0.0.1a',
    description='Neat molecular geometry tools',
    author='Niccolo Ricardi',
    author_email='Niccolo.Ricardi@unige.ch',
    package_dir={'geomtools': 'geomtools'},
    packages=['geomtools'],
)