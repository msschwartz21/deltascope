from setuptools import setup

def readme():
	with open('README.rst') as f:
		return (f.read())

def requirements():
	with open('requirements.txt') as f:
		return (f.read())

setup(name='cranium',
	version='0.1.4',
	description='Python package to quantify biological structure',
	long_description=readme(),
	keywords='biology image analysis quantification',
	url='https://github.com/msschwartz21/craniumPy',
	author='Morgan Schwartz',
	author_email='mschwartz@smith.edu',
	license='GNU',
	packages=['cranium'],
	install_requires=[
		'cycler==0.10.0',
		'decorator==4.1.2',
		'h5py==2.7.1',
		'matplotlib==1.5.0',
		'networkx==2.0',
		'numpy==1.14.0',
		'pandas==0.22.0',
		'patsy==0.4.1',
		'Pillow==5.0.0',
		'pyparsing==2.2.0',
		'python-dateutil==2.6.1',
		'pytz==2017.3',
		'PyWavelets==0.5.2',
		'scikit-image==0.13.1',
		'scikit-learn==0.19.1',
		'scipy==1.0.0',
		'six==1.11.0',
		'statsmodels==0.8.0'
	],
	include_package_data=True
	)