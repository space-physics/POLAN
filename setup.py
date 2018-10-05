#!/usr/bin/env python
import setuptools  # noqa: F401
from numpy.distutils.core import setup, Extension
from pathlib import Path
import os


if os.name == 'nt':
    sfn = Path(__file__).parent / 'setup.cfg'
    stxt = sfn.read_text()
    if '[build_ext]' not in stxt:
        with sfn.open('a') as f:
            f.write("[build_ext]\ncompiler = mingw32")

setup(ext_modules=[Extension(name='runpolan',
                             sources=['src/polan.f', 'src/polmis.f',
                                      'src/polsin.f', 'src/polsub.f', 'src/polrunsub.f'],
                             extra_f77_compile_args=['-Wno-line-truncation'])],
      )
