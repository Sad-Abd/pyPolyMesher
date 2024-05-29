from setuptools import setup, find_packages

with open("README.md", "r", encoding = "utf-8") as fh:
    long_description = fh.read()

setup(
    name="pyPolyMesher",
    version="1.0.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    author="Sad-Abd",
    author_email="abedisadjad@gmail.com",
    description="A Python package for polygonal mesh generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Sad-Abd/pyPolyMesher",
    license="GPL-3.0",
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "ezdxf"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11"
    ],
)