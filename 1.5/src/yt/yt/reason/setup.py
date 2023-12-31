#!/usr/bin/env python
import setuptools
import os, sys, os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('reason',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    return config
