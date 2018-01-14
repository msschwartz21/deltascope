from setuptools import setup

def readme():
	with open('README.rst') as f:
		return (f.read())

def requirements():
	with open('requirements.txt') as f:
		return (f.read())

setup(name='cranium',
	version='0.2.1',
	description='Python package to quantify biological structure',
	long_description=readme(),
	keywords='biology image analysis quantification',
	url='https://github.com/msschwartz21/craniumPy',
	author='Morgan Schwartz',
	author_email='mschwartz@smith.edu',
	license='GNU GPL',
	packages=['cranium'],
	install_requires=[
		'cycler==0.10',
		'decorator==4.1',
		'h5py==2.7',
		'matplotlib>=2.0,<3.0',
		'networkx==2.0',
		'numpy>=1.0,<2.0',
		'pandas==0.22',
		'patsy==0.4',
		'Pillow==5.0',
		'pyparsing==2.2',
		'python-dateutil==2.6',
		'pytz==2017.3',
		'PyWavelets==0.5',
		'scikit-image==0.13',
		'scikit-learn==0.19',
		'scipy==1.0',
		'six==1.11',
		'statsmodels==0.8'
	],
	include_package_data=True
	)