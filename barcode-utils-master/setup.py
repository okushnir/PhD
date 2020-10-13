import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Barcode-utils", # Replace with your own username
    version="0.0.1",
    author="Maoz Gelbart",
    author_email="author@example.com",
    description="A package with utilities to handle reads with molecular identifiers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MaozGelbart/mypackage",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)