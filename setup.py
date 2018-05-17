#!/usr/bin/env python
install_requires=['numpy','xarray']
tests_require=['pytest','nose','coveralls']
# %%
from setuptools import find_packages
from numpy.distutils.core import setup,Extension

setup(name='polan',
      packages=find_packages(),
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
            sources=['polan.f', 'polmis.f', 'polsin.f', 'polsub.f','polrunsub.f'],
                    f2py_options=['--quiet'],
                    extra_f77_compile_args=['-Wno-line-truncation'])],
       install_requires=install_requires,
       tests_require=tests_require,
       extras_require={'tests':tests_require},
       python_requires='>=3.5',
	  )
