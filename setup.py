import numpy as np
from Cython.Build import cythonize
from setuptools import setup, Extension

extensions = cythonize(
    [
        Extension("sequence", ["src/fragments_h5/sequence.pyx"]),
    ],
    compiler_directives={"language_level": "3"},
)

setup(
    include_dirs=[np.get_include()],
    ext_modules=extensions,
)
