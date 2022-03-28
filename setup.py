import setuptools

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scanexitronlr",
    version="1.1.7",
    author="Josh Fry",
    author_email="fryxx094@umn.edu",
    description="ScanExitronLR: a lightweight tool for the "
                "characterization and quantification of exitrons "
                "in long read RNA-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ylab/ScanExitronLR",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.7",
    ],
    install_requires=["pysam", "liqa", "gffutils", "pandas", "biopython"],
    python_requires='>=3.7',
    packages=["src"],
    package_dir={"src": "src"},
    # include_package_data=True,
    package_data={'': ['*.tsv']},
    entry_points={
        "console_scripts": [
            "selr=src.__main__:main",
            "scanexitronlr=src.__main__:main"
        ]
    },
)
