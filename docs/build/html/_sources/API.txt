.. _api:
API
====

.. class:: brain(filepath)

	Object to manage biological data and associated functions. Initialization of this object automatically reads the hdf5 file and selects the appropriate channel from the Probability file and saves as global :py:attr:`brain.raw_data` attribute.

	:param str filepath: Filepath to hdf5 probability file
	:return: :py:attr:`brain.raw_data`
	:rtype: 3D array

	.. py:attribute:: brain.raw_data

		Array of shape [z,y,x] containing raw probability data

	.. py:attribute:: brain.df

		Dataframe with four columns: x,y,z,value with all points in :py:attr:`brain.raw_data`

	.. py:attribute:: brain.threshold

		Value used to threshold the data prior to calculating the model

	.. py:attribute:: brain.df_thresh

		Data frame containing only points with values above the specified threshold

	.. py:attribute:: brain.model

		Dictionary containing the coefficients for the parabolic model of the data

.. py:method:: brain.read_data(filepath)

	Reads 3D data from file and selects appropriate channel based on the assumption that the channel with the most zeros has zero as the value for no signal

	:param str filepath: Filepath to hdf5 probability file
	:return: Array of shape [z,y,x] containing raw probability data

.. py:method:: brain.show_plane(dimension,plane)

	Displaces single specified plane in a matplotlib figure

	:param str dimension: string specifying dimension of plane, e.g. x,y,z
	:param int plane: Integer specifying position of plane
	:return: Matplotlib figure

.. py:method:: brain.create_dataframe()

	Creates a pandas dataframe containing the x,y,z and signal/probability value for each point in the :py:attr:`brain.raw_data` array

	:return: :py:attr:`brain.df`
	:rtype: Dataframe with four columns: x,y,z,value

.. py:method:: brain.fit_model(threshold)

	Calculates the mathematical model of the data by identifying the flat plane and parabolic plane that can fit the data before calculating their intersect

	:param float threshold: float value between 0 and 1, used to select lower bound of values
	:return: :py:attr:`brain.model`
	:rtype: dictionary containing the coefficients for the parabolic model of the data