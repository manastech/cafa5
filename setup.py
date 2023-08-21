from setuptools import setup, find_packages

setup(
    name="manas_cafa",
    version="0.0.0",
    packages=["manas_cafa"],
    python_requires=">=3.8",
    package_data={
        '': ['LICENSE', '*.md'],
    },
    install_requires=[
        'bio >=1.5.9, <2.0.0',
        'numpy >=1.22, <1.24',
        'pandas >=2.0.1, <3.0.0',
        'tensorflow >=2.8.2, <3.0.0',
        'networkx >=3.1, <4.0',
        'obonet >=1.0.0, <2.0.0',
    ],
    description="manas ML model for cafa5 contest",
    url="https://github.com/manastech/cafa5",
)
