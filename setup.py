from setuptools import setup

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='cranium',
	version='0.0.1',
	description='Python package to quantify biological structure',
	long_description=readme(),
	classifiers=[
		'License :: MIT License'
		'Programming Language :: Python :: 3.5.2',
		'Topic :: Biology :: Image analysis'
		],
	keywords='biology image analysis quantification',
	url='https://github.com/msschwartz21/craniumPy',
	author='Morgan Schwartz',
	author_email='mschwartz@smith.edu',
	license='MIT',
	packages=['cranium'],
	install_requires=[
		'os'
		],
	include_package_data=True
	)