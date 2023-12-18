from setuptools import setup

setup(
    name="fragment_h5",
    version="1.0",
    packages=['fragments_h5'],
    scripts=["bin/build-fragments-h5"],
    install_requires=[
        "numpy",
        "h5py",
        "pysam",
        "funcsigs",
    ]
)
