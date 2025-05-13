from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pyMut",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Una librería Python para visualizar mutaciones genéticas",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/pyMut",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "numpy>=1.20.0",
        "seaborn>=0.11.0",
    ],
    include_package_data=True,
    package_data={
        "pyMut": ["data/examples/*.tsv"],
    },
) 