from setuptools import setup

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

setup(
    name = "sign_vectors",
    version = readfile("VERSION").strip(), # the VERSION file is shared with the documentation
    description='A package for SageMath to work with sign vectors and oriented matroids',
    long_description = readfile("README.md"),
    long_description_content_type="text/markdown",
    url='https://github.com/MarcusAichmayr/sign_vectors',
    author='Marcus Aichmayr',
#    author_email='',
    license='GPLv3',
    classifiers=[
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 5 - Stable',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: GNU General Public License v3',
      'Programming Language :: Python :: 3.8.5',
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords = ['sign vectors', 'elementary vectors', 'oriented matroids'],
    packages = [
        'elementary_vectors',
        'sign_vectors',
        'sign_vectors.oriented_matroids',
    ],
#    setup_requires   = ['sage-package'],
#    install_requires = ['sage-package', 'sphinx'],
)

