from setuptools import setup

setup(name='pydebruijn',
      version='0.1',
      description='A library for working on de Bruijn sequences and related concepts.',
      url='http://github.com/adamasstokhorst/pydebruijn',
      author='Adamas Aqsa Fahreza',
      author_email='adamas@ntu.edu.sg',
      license='GPL',
      packages=['pydebruijn'],
      install_requires=['sympy', 'networkx'],
      include_package_data=True,
      zip_safe=False)
