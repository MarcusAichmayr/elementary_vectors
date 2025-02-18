from setuptools import setup


# Get information from separate files (README)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


setup(
    name="elementary_vectors",
    version="v1.2",
    description="a SageMath package to work with elementary vectors, sign vectors, oriented matroids and vectors with components in intervals",
    long_description=readfile("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/MarcusAichmayr/elementary_vectors",
    author="Marcus S. Aichmayr",
    author_email="aichmayr@mathematik.uni-kassel.de",
    license="GPLv3",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords=[
        "elementary_vectors",
        "vectors",
        "intervals",
        "linear inequality systems",
        "sign vectors",
        "oriented matroids",
    ],
    packages=[
        "elementary_vectors",
        "vectors_in_intervals",
        "sign_vectors",
    ],
    extras_require={
        "passagemath": [
            "passagemath-symbolics",
            "passagemath-flint",
            "passagemath-graphs",
            "passagemath-repl",
        ],
    },
)
