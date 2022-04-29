import setuptools
import re
import os
import codecs

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

PKG = "spcount"

version = find_version("src/spcount/__version__.py")

setuptools.setup(
    name=PKG,
    version=version,
    author="Quanhu Sheng",
    author_email="quanhu.sheng.1@vumc.org",
    description="Multiple genomes read count",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shengqh/spcount",
    download_url="https://github.com/shengqh/spcount/archive/v" + version + ".tar.gz",
    project_urls={
        "Bug Tracker": "https://github.com/shengqh/spcount/issues",
    },    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src", exclude=["tests", "tests.*", "scripts"]),
    python_requires='>=3.6',
    entry_points = {
        'console_scripts': ['spcount=spcount.__main__:main'],
    },
    package_data={'': ['src/spcount/slurm.template']},
    install_requires=['argparse', 'pandas', 'seaborn' ],
    include_package_data=True,
    zip_safe=False,
)

