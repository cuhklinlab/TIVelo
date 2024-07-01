from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="tivelo",
    version="0.0.3",
    description="single cell velocity analysis",
    package_dir={"": "tivelo"},
    packages=find_packages(where="tivelo"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aqlkzf/tivelo",
    author="Muyang GE",
    author_email="muyang@mail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    install_requires=["bson >= 0.5.10"],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2"],
    },
    python_requires=">=3.9",
)