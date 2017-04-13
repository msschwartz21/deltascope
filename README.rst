craniumPy
##########

This package is designed to quantify biological structures in 3D image data. 

Installation
++++++++++++

Package hosted on `PyPI <https://pypi.python.org/pypi/cranium>`_. ::

	$ pip install cranium

Documentation
++++++++++++++

Complete documentation is available on `Read the Docs <http://craniumpy.readthedocs.io/en/latest/>`_.

Basic Usage
++++++++++++

.. code-block:: python

	import(cranium)
	file = "C:\\Users\\zfishlab\\Desktop\\zrf1wt13umyot21um\\wt\\AT2\\AT_04_Probabilities.h5"

	sample = cranium.brain(file)
	sample.create_dataframe()
	sample.fit_model(0.9)
	sample.plot_model()