from setuptools import setup


def readme():
    with open('README.md', 'rt') as f:
        return f.read()


setup(name='ngmaster',
      version='0.5',
      description='In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GPLv2',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics Neisseria sequence typing',
      url='https://github.com/MDU-PHL/ngmaster',
      author='Jason Kwong',
      author_email='kwongj@gmail.com',
      license='GPLv2',
      packages=['ngmaster'],
      install_requires=[
          'argparse',
          'BioPython',
          'bs4',
          'requests',
      ],
      test_suite='nose.collector',
      tests_require=[],
      entry_points={
          'console_scripts': ['ngmaster=ngmaster.ngmaster:main'],
      },
      include_package_data=True,
      zip_safe=False)
