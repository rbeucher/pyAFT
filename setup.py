from setuptools import setup, find_packages, Extension


module1 = Extension('ketcham',
                    include_dirs=['include'],
                    sources = ['./src/Ketcham.c'])


setup(name='PyAFT',
      version='0.1',
      description='Apatite Fission Track utilities',
      url='',
      author='Romain Beucher',
      author_email='romain.beucher@geo.uib.no',
      license='MIT',
      packages=find_packages(),
      install_requires=['numpy','matplotlib'],
      zip_safe=False,
      ext_modules=[module1]
      )

