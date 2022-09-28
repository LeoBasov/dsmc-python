from setuptools import setup
setup(
    name='dsmc',
    version='0.6.0',
    author='Leo Basov',
    python_requires='>=3.6, <4',
    packages=["dsmc"],
    install_requires=[
        'numpy',
        'llvmlite',
        'scipy',
        'numba'
    ],
)
