from setuptools import setup
setup(
    name='iq'
    version='0.0.1',
    author='',
    maintainer='',
    description='A python package for protein quantification',
    long_description="This R package provides an implementation of the MaxLFQ algorithm by Cox et al. (2014) in a comprehensive pipeline for DIA-MS (Pham et al. 2020). It also offers options for protein quantification using the N most intense fragment ions, using all fragment ions, and the Tukey's median polish algorithm. In general, the tool can be used to integrate multiple proportional observations into a single quantitative value.",
    url='',
    download_url='',
    keywords='preprocess, create_protein_list, maxLFQ, create_protein_table',
    python_requires='>=3.7, <4',
    packages=find_packages(include=['maxLFQ', 'maxLFQ.*']),
    install_requires=[
        'pandas==0.23.3',
        'numpy>=1.14.5',
        'matplotlib>=2.2.0'
    ],
    package_data={},
    entry_points={}
)