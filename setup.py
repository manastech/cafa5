from setuptools import setup, find_packages

setup(
    name="manas-cafa5",
    version="0.0.0",
    packages=[],
    python_requires=">=3.9",
    package_data={
        '': ['LICENSE', '*.md'],
    },
    install_requires=[
        'bio >=1.5.9, <2.0.0',
        'numpy >=1.22, <1.24',
        'pandas >=2.0.1, <3.0.0',
        'tensorflow >=2.12.0, <3.0.0',
    ],
    description="manas ML model for cafa5 contest",
    url="https://github.com/manastech/cafa5",
)
