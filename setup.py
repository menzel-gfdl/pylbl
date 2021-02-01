from setuptools import find_packages, setup


with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setup(
    name="pyrad",
    version="0.0.1",
    author="R. Menzel",
    author_email="author@example.com",
    description="A is a simple (one-dimensional) all-sky atmospheric radiation package.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/menzel-gfdl/pylbl",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: LGPL-3.0 License",
        "Operating System :: OS Independent",
    ],
    install_requires=["netCDF4", "numpy", "scipy"],
    python_requires=">=3.5")
