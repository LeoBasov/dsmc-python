from setuptools import setup
setup(
    name='dsmc',
    version='0.9.0',
    author='Leo Basov',
    python_requires='>=3.10, <4',
    packages=["dsmc", "dsmc.writer", "dsmc.mesh"],
    install_requires=[
        'numpy',
        'llvmlite',
        'scipy',
        'numba'
    ],
)
