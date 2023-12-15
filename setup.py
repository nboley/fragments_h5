from setuptools import setup, find_packages

setup(
    name="fragment_h5",
    version="1.0",
    packages=find_packages(),
    scripts=["bin/build_fragment_h5.py"],
)
