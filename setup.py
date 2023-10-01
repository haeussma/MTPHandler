from setuptools import setup, find_packages

setup(
    name="MTPHandler",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["numpy", "sympy", "pandas", "matplotlib", "plotly"],
    extras_require={
        "JupyterNotebook": ["jupyter", "ipykernel"],
    },
)
