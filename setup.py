from setuptools import setup, find_packages

setup(
    name='mtphandler',
    version='0.1.0',
    packages=find_packages(include=['PlateHandler', 'PlateHandler.*']),
    install_requires=[
        'numpy',
        'sympy',
        'pandas',
        'matplotlib',
        'plotly'
    ],
    extras_require={
        'JupyterNotebook': ['jupyter', 'ipykernel'],
    }


)
