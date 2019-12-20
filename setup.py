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
    packages=["bmeg", "bmeg.util", "bmeg.enrichers"],
    package_data={"": ["bmeg-dictionary/gdcdictionary/schemas/*.yaml"]},
    zip_safe=False,
    python_requires=">=3.6, <4",
    install_requires=[
        "dataclasses>=0.6",
        "jsonschema>=3.0.1",
        "requests>=2.19.1",
        "requests-cache>=0.4.13",
        "dictionaryutils==2.0.7"
    ],
    extras_require={
        "test": [
            "nose>=1.3.7",
            "flake8>=3.5.0",
        ]
    },
    dependency_links=[
        "git+https://github.com/uc-cdis/dictionaryutils.git@2.0.7"
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: MIT License",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: Python :: 3.7",
    ],
)
