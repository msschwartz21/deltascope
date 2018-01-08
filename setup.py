from setuptools import setup

def readme():
	with open('README.rst') as f:
		return f.read()

def requirements():
	with open('requirements.txt') as f:
		return f.read()

setup(name='cranium',
	version='0.1.2',
	description='Python package to quantify biological structure',
	long_description=readme(),
	keywords='biology image analysis quantification',
	url='https://github.com/msschwartz21/craniumPy',
	author='Morgan Schwartz',
	author_email='mschwartz@smith.edu',
	license='GNU',
	packages=['cranium'],
	install_requires=requirements(),
	include_package_data=True
	)