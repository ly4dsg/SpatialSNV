from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="spatialsnv",  
    version="1.0.2", 
    description="A toolkit for spatial SNV analysis", 
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Yi Liu",
    author_email="liuyi6@genomics.cn", 
    url="https://https://github.com/ly4dsg/SpatialSNV",  
    license="MIT", 
    packages=find_packages(), 
    py_modules=["spatialSNV"],
    python_requires=">=3.7",  
    install_requires=[
        "multiprocess==0.70.17",
        "click==8.1.8",
        "pysam==0.22.1",
        "pandas==1.5.3",
        "numpy==1.24.4",
        "scikit-learn==1.4.1.post1",
        "scipy==1.11.4",
        "igraph==0.11.5",
        "seaborn==0.13.2",
        "matplotlib==3.7.5",
        "scanpy==1.10.0",
        "tqdm",
        "cairocffi",
        "pycairo",
        "leidenalg",
    ],
    entry_points={
        "console_scripts": [
            "spatialsnvtools=spatialsnvtools.cli:main",  
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
