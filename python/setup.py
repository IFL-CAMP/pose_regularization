from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np
import os
from io import open  # for python 2

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, '../README.md'), encoding='utf-8') as f:
    long_description = f.read()
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description += f.read()

setup(
    name='pose_regularization',
    license='LGPLv3',
    version='0.1.0',
    description='Total Variation Regularization of Pose Signals',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/IFL-CAMP/pose_regularization',
    author='Marco Esposito',
    author_email='marco.esposito@tum.de',
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(
        "pose_regularization/pypose_regularization.pyx",
        include_path=['../include', np.get_include()]
    ),
    packages=['pose_regularization'],
    install_requires=[
        'enum34;python_version<"3.4"',
        'cython',
        'numpy',
    ],
    tests_require=[
        'cython',
        'numpy',
        'scipy',
    ]
)
