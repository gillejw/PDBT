from setuptools import setup, find_packages

# read the contents of the README.md file
from os import path
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pdbt',
    use_scm_version=True,
    description='PDBT - Phage Display Bioinformatics Toolkit',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='http://github.com/gillejw/PDBT',
    project_urls={
        "Issues": "https://github.com/gillejw/PDBT/issues"
    },
    author='James W. Gillespie',
    author_email='gillejw@auburn.edu',
    license='MIT',
    keywords=['phage display'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages=find_packages(include=['pdbt']),
    install_requires=[],
    extras_require={
        'dev': []
    },
    test_suite='pytest',
    tests_require=['pytest'],
    setup_requires=['flake8','setuptools_scm'],
    include_package_data=True,
    python_requires='>=3',
    zip_safe=False
    )
