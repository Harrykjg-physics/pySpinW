try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

setup(
    name='pySpinW',
    version='0.1',
    packages=find_packages(),
    url='http://www.spinw.org',
    license='MIT License',
    author='Simon Ward',
    author_email='simon.ward@esss.se',
    description="Python implementation of SpinW",
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
    ],
    include_package_data=True,
    install_requires=['numpy']
)
