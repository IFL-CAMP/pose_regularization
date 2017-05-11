from setuptools import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    name='pose_denoising',
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("pose_denoising",
                  sources=["pose_denoising.pyx", "../src/manifoldGeneric.cpp", "../src/manifoldSE3prod.cpp"],
                  include_dirs=[numpy.get_include(), '../include'],
                  language="c++",  # generate C++ code
                  extra_compile_args=["-std=c++11", "-DSE3PROD"],
                  extra_link_args=["-std=c++11"]
                  )
    ]
)
