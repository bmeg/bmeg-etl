import io
import os

from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


setup(
    name="bmeg",
    version="0.1.0",
    description="Models and utilities for BMEG-ETL processes",
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url="https://github.com/bmeg/bmeg-etl",
    license="MIT",
    package_dir={"": "src/"},
    packages=["bmeg"],
    python_requires=">=3.6, <4",
    install_requires=[
        "dataclasses>=0.6",
        "mygene>=3.0.0",
    ],
    extras_require={
        "test": [
            "nose>=1.3.7",
            "flake8>=3.5.0",
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: MIT License",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)
