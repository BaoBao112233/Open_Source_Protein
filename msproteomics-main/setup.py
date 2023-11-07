from setuptools import setup

with open('README.md') as f:
    long_description = f.read()

setup(
    name='msproteomics',
    version='0.0.2',    
    description='A Python package for mass spectrometry-based proteomics data processing',
    url='https://github.com/tvpham/msproteomics',
    author='Thang Pham',
    author_email='t.pham@amsterdamumc.nl',
    license='Apache License 2.0',
    packages=['msproteomics'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
		'License :: OSI Approved :: Apache Software License',
    ],
)
