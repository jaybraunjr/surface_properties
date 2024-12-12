from setuptools import setup, find_packages

setup(
    name='surface_properties',  # Replace with your package name
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'MDAnalysis',
        'matplotlib',
    ],
    author='Your Name',
    author_email='u1266568@utah.edu',
    description='A package for molecular dynamics analysis',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
