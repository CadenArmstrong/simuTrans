#! /usr/bin/env python

from distutils.core import *
from distutils import sysconfig
from distutils.extension import Extension

import numpy

numpy_include = numpy.get_include()
_transitmodel = Extension("_transitmodel",
        ["transitmodel_wrap.cxx",
            "transitmodel.cpp"],
        include_dirs = [numpy_include],
        )

setup(name="Simutrans",
        author_email = "",
        url = "",
        version = "0.0.0",
        py_modules = ["transitmodel"],
        ext_modules = [_transitmodel])

