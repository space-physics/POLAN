#!/usr/bin/env python
req=['python-dateutil','pytz','nose','numpy','xarray','matplotlib','seaborn']
# %%
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception as e:
    import pip
    pip.main(['install'] + req)
# %%
import setuptools #enables develop
from numpy.distutils.core import setup,Extension

setup(name='polan',
      packages=['polan'],
      author='Michael Hirsch, Ph.D',
      description='Model of Earth ionosphere true height.',
      version='0.1.0',
      url='https://github.com/scivision/polan',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 3 - Alpha',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
      ext_modules=[Extension(name='runpolan',
            sources=['polan.f', 'polmis.f', 'polsin.f', 'polsub.f'],
                    f2py_options=['--quiet'],
                    extra_f77_compile_args=['-Wno-line-truncation'])]
	  )
