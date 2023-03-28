from setuptools import setup, find_packages

classifiers = [
  'Development Status :: 5 - Production/Stable',
  'Intended Audience :: Education',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]

setup(
  name='crystal22',
  version='0.0.1',
  description='A library for running and plotting crystal22 structures',
  long_description=open('README.txt').read() + '\n\n' + open('CHANGELOG.txt').read(),
  url='',
  author='William Comaskey',
  author_email='williamcomaskey@gmail.com',
  license='MIT',
  classifiers=classifiers,
  keywords='crystal22',
  packages=find_packages(),
  install_requires=['']
)
