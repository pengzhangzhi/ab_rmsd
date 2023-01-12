from setuptools import setup, find_packages

__version__ = "0.0.1"
setup(
  name = 'antibody-RMSD',
  packages = find_packages(),
  version = __version__,
  license='MIT',
  description = 'Calculate RMSD between two antibody structures',
  author = 'Zhangzhi Peng',
  author_email = 'pengzhangzhics@gmail.com',
  url = '',
  long_description_content_type = 'text/markdown',
  keywords = [
    'protein',
    'antibody'
  ],
  install_requires=[
    "torch>=1.7",
    "easydict",
    "biopython",
  ],
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.6',
  ],
  scripts=['abrmsd'],
)