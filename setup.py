#Copyright (c) 2018, BGI-Shenzhen

from setuptools import setup

def readme():
    with open("README.rst") as f:
        return f.read()

setup(name='EPIC',
      version='1.0',
      description='Neoantigen prediction tool',
      long_description=readme(),
      entry_points = {'console_scripts': ['epic-predict=epic.command_line:main']},
      url='https://gitlab.genomics.cn/cell_therapy/EPIC',
      author='Weipeng Hu',
      author_email='huweipeng@genomics.cn',
      license='MIT',
      packages=['epic'],
      install_requires=[
        'argparse',
        'numpy',
        'pandas',
        'scikit-learn',
        'setuptools',
      ],
      include_package_data=True,
      zip_safe=False)
