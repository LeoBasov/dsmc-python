from setuptools import setup
setup(
    name='dsmc',
    version='0.7.0',
    author='Leo Basov',
    python_requires='>=3.10, <4',
    packages=["dsmc"],
    install_requires=[
        'numpy',
        'llvmlite',
        'scipy',
        'numba'
    ],
)
