import os
import itertools
from setuptools import setup, find_packages


def get_long_description():
    """Load the contents of README.md to use as a long description
    
    Returns:
        The content of README.md as a string.
    """
    this_directory = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
        long_description = f.read()
    return long_description


def main():
    # Package data
    PACKAGE_DATA = {
        "name": "ManipulateAggregates",
        "version": "0.1.3",
        "description": "Manipulate molecular DOF and scan PES of aggregates",
        "author": "Torsten Sachse",
        "mail": "torsten.sachse@gmail.com",
        "url": "https://github.com/razziel89/ManipulateAggregates",
        "long_description": get_long_description(),
        "long_description_content_type": "text/markdown",
        "classifiers": [
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
            "Natural Language :: English",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: C++",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3 :: Only",
        ],
        "install_requires": [
            "numpy >=1.10",
            "FireDeamon >=0.1.2,<0.2",
            "maagbel >=0.1.2,<0.2",
        ],
        "package_data": {"": ["data/*"]},
    }

    setup(
        packages=find_packages(),
        entry_points={
            "console_scripts": [
                "energyscan = ManipulateAggregates.bin.energyscan:entrypoint",
                "hashsort = ManipulateAggregates.bin.hashsort:entrypoint",
                "manipagg = ManipulateAggregates.bin.manipagg:entrypoint",
            ]
        },
        **PACKAGE_DATA,
    )


if __name__ == "__main__":
    main()
