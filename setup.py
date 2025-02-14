from setuptools import setup, find_packages

setup(
    name="sprint-genome-analyzer",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "networkx>=2.6.0",
        "scipy>=1.7.0",
        "biopython>=1.79",
        "requests>=2.26.0",
        "dask>=2021.8.0",
        "pyarrow>=5.0.0",
        "scikit-learn>=0.24.0",
        "torch>=1.9.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.5b2",
            "flake8>=3.9.0",
            "mypy>=0.910",
        ]
    },
    python_requires=">=3.8",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive genome analyzer for sprint performance",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    keywords="genomics, bioinformatics, sports-science, network-analysis",
)
