from setuptools import setup, find_packages

setup(
    name='Monte_Carlo_Shielding',
    version='0.1',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)
