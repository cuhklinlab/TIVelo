from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="tivelo",
    version="0.0.5",
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
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent",
    ],
    install_requires=["bson >= 0.5.10",
                        "numpy >= 1.21.0",
                        "scipy >= 1.10.0", 
                         "matplotlib >= 3.4.0",
                        "pandas >= 2.0.0",
                        "networkx >= 3.0.0",
                        "scikit-learn >= 1.0.0",
                        "scanpy >= 1.9.0",
                        "tqdm >= 4.0.0",
                        "torch >= 1.9.0",
                        "scvelo >= 0.3.0",
                        "numba >= 0.55.0",
                        
                        
                          
                      ],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2"],
    },
    python_requires=">=3.9",
)