from setuptools import setup, find_packages

setup(
    name="tivelo",
    version="0.0.7",
    packages=find_packages(include=['tivelo', 'tivelo.*']),
    include_package_data=True,
)
