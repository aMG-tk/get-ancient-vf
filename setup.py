from setuptools import setup
import versioneer

requirements = [
    "pandas>=1.3.4",
    "scipy>=1.7.3",
    "tqdm>=4.62.3",
    "numpy>=1.21.4",
    "biopython>=1.79",
    "requests>=2.26.0",
]

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        "setuptools>=18.0",
        "Cython>=0.29.24",
    ],
    name="get-ancient-vf",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple tool to map ancient reads against VFDB.",
    license="GNUv3",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/get-ancient-vf",
    packages=["get_vf"],
    entry_points={"console_scripts": ["getVF=get_vf.__main__:main"]},
    install_requires=requirements,
    keywords="get-ancient-vf",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
