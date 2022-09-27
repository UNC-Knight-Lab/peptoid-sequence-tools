from setuptools import setup, find_packages

VERSION = '1.0.1'
DESCRIPTION = 'My first Python package'
LONG_DESCRIPTION = 'My first Python package with a slightly longer description'

# Setting up
setup(
    name="seq_match",
    version=VERSION,
    author="Supraja Chittari, Alan Wang",
    author_email="suprajac@email.unc.edu, alanwang205@gmail.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url='https://github.com/UNC-Knight-Lab/peptoid-sequence-tools',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'xlsxwriter'
    ],
    entry_points={
        'console_scripts': ['seq_match=seq_match.command_line:main'],
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
    ],
    license='GPL-3.0'
)