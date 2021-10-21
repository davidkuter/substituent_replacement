from setuptools import setup
from setuptools import find_packages

# App dependencies to for being imported into other applications.
app_requirements = [
    'loguru==0.5.*',
    'xmltodict==0.12.*',
]

# Dev dependencies only required for development
dev_requirements = [
    'pytest==5.4.*',
    'pytest-cov==2.8.*',
    'pytest-mock==3.1.*',
    'pytest-html>=3.1.*',
    'pylint==2.6.*',
    'ipython==7.14.*',
    'tox==3.15.*',
    'pip-tools==6.4.*',
]

# extract version tag from code
version = dict()
with open('src/substituent_replacement/__init__.py') as fn:
    exec(fn.read(), version)

setup(name='substituent_replacement',
      description='Substituent Replacement API',
      long_description='Method to generate derivatives of a molecule based on frequently used R-group transformations',
      version=version['__version__'],
      url='https://github.com/davidkuter/substituent_replacement',
      license='CC-BY',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      zip_safe=False,
      install_requires=app_requirements,
      entry_points={},
      extras_require={'dev': dev_requirements})
