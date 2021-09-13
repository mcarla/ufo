import setuptools


setuptools.setup(
    name="ufo",
    version="0.01",
    author="Michele Carla'",
    author_email="mcarla@cells.es",
    description="Unreliable but Fast Optics code",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="BSD",
    url="https://github.com/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)

