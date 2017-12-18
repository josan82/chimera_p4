#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
# Python stdlib
import os
from contextlib import contextmanager


@contextmanager
def ignored(*exceptions):
    """
    Ignore all exceptions
    """
    try:
        yield
    except exceptions:
        pass


@contextmanager
def enter_directory(path):
    """
    Change directory to `path` and then go back to original one
    """
    oldwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(oldwd)