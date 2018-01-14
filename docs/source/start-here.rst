.. _start here: 

Start Here
===========

.. _right for you:

Is Cranium right for you? 
+++++++++++++++++++++++++++

Cranium may be able to help you if:

- your data consists of a set of 3D image stacks
- your data contains a clear structure and shape that has consistent gross morphology between control and experimental samples
- you want to identify extreme or subtle differences in the structure between your experimental groups

Cranium cannot help if:

- your data was collected using cryosections that need to be aligned after imaging
- your experiment changes the gross anatomy of the structure between control and experimental samples

.. Todo:: Data format requirements
	-ImageJ created H5 images
	-No more than 4 different image channels 
	-Ect

.. _install:

Installation
+++++++++++++

Cranium is most easily installed using `pip`_, Python's package installer. However, Cranium's dependencies can be difficult to install using pip, so we recommend that you install `Anaconda`_ to manage more complicated packages. See `Setting up a Python environment <python set up>`_ for more details. 

.. Todo:: Cranium can be installed using `pip`_, Python's package installer. Some Cranium dependencies are not available through pip, so we recommend that you install `Anaconda`_ which automatically includes and installs the remaining dependencies. See `Setting up a Python environment <python set up>`_ for more details.

.. code::
	
	$ pip install cranium

.. note:: If you are unfamiliar with the command line, check out `this tutorial <command line tutorial>`_. Additional resources regarding the command line are available `here <resources>`_.

.. warning:: Packages required for cranium depend on Visual C++ Build Tools, which can be downloaded `here <Build Tools>`_.

.. _python set up:

Setting up a Python environment
++++++++++++++++++++++++++++++++

If you're new to scientific computing with Python, we recommend that you install `Anaconda`_ to manage your Python installation. Anaconda is a framework for scientific computing with Python that will install important packages (`numpy`_, `scipy`_, and `matplotlib`_).

.. warning:: Cranium is written in Python 3, and requires the installation the Python 3 version of Anaconda. 

.. _Anaconda: https://www.anaconda.com/what-is-anaconda/
.. _pip: https://en.wikipedia.org/wiki/Pip_(package_manager)
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _matplotlib: https://matplotlib.org/
.. _command line tutorial: http://www.vikingcodeschool.com/web-development-basics/a-command-line-crash-course
.. _pip tutorial: https://programminghistorian.org/lessons/installing-python-modules-pip
.. _Build Tools: http://landinghub.visualstudio.com/visual-cpp-build-tools
