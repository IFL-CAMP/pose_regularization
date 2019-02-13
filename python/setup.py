from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np

setup(
    name='pose_regularization',
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
