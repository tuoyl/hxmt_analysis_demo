import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hxmt_analysis_demo",
    version="0.1.0",
    author="Youli Tuo",
    author_email="tuoyl@ihep.ac.cn",
    description="a package for hxmt analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tuoyl/hxmt_analysis_demo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
