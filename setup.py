from setuptools import setup, find_packages

setup(
    name="biomatch",
    version="0.5.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pyfaidx",
        "pandas",
        "numpy",
        "biopython",
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            "biomatch=biomatch.biomatch:main",
        ],
    },
    package_data={
        "biomatch": ["analysis_scripts/*", "kmer_ref_panels/*"],
    },
    python_requires=">=3.6",
    description="BioMatch: A data-driven sample identification matcher",
    author="VonPoo",
    author_email="fengbobo927@163.com",
    url="https://github.com/vonpoo/BioMatch",
)