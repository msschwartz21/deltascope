from setuptools import setup

def readme():
	with open('README.rst') as f:
		return (f.read())

def requirements():
	with open('requirements.txt') as f:
		return (f.read())

setup(name='deltascope',
	version='0.2.3',
	description='Python package to quantify biological structure',
	long_description=readme(),
	keywords='biology image analysis quantification',
	url='https://github.com/msschwartz21/deltascope',
	author='Morgan Schwartz',
	author_email='msschwartz21@gmail.com',
	license='GNU GPL',
	packages=['deltascope'],
	install_requires=[
		'cycler==0.10',
		'decorator==4.1',
		'h5py==2.8',
		'matplotlib>=1.5,<3.0',
		'networkx==1.11',
		'numpy>=1.0,<2.0',
		'pandas>=0.19',
		'patsy==0.4',
		'Pillow==5.0',
		'pyparsing==2.2',
		'python-dateutil==2.6',
		'pytz==2017.3',
		'PyWavelets==0.5',
		'scikit-image==0.12',
		'scikit-learn==0.18',
		'scipy>=0.18',
		'six==1.11',
		'statsmodels==0.8'
	],
	include_package_data=True
	)
