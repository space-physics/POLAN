#!/usr/bin/env python
import setuptools  # noqa: F401
from numpy.distutils.core import setup, Extension

setup(
    ext_modules=[
        Extension(
            name="runpolan",
            sources=[
                "src/polan.f",
                "src/polmis.f",
                "src/polsin.f",
                "src/polsub.f",
                "src/polrunsub.f",
            ],
            extra_f77_compile_args=["-std=legacy"],
        )
    ],
)
