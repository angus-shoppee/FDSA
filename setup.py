# Feature-Directed Splice Analysis
# Copyright (C) 2024 Angus Shoppee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from setuptools import setup, find_packages


setup(
    name="fdsa",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    py_modules=["main"],
    install_requires=[
        "conorm>=1.2.0",
        "entrezpy>=2.1.3",
        "gtfparse>=1.2.1",
        "matplotlib>=3.0.0",
        "numpy>=1.22.2",
        "pandas>=1.2.4",
        "pybiomart>=0.2.0",
        "pysam>=0.15.4"
    ],
    entry_points={
        "console_scripts": [
            "fdsa=main:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ]
)
