.. _api:
API
====

.. class:: brain(filepath)

	Object to manage biological data and associated functions. 

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

	.. py:attribute:: brain.subset

		Random sample of the dataframe from :py:attr:`brain.df_thresh`

.. py:method:: brain.read_data(filepath)

	Reads 3D data from file and selects appropriate channel based on the assumption that the channel with the most zeros has zero as the value for no signal

	:param str filepath: Filepath to hdf5 probability file
	:return: Array of shape [z,y,x] containing raw probability data

.. py:method:: brain.create_dataframe()

	Creates a pandas dataframe containing the x,y,z and signal/probability value for each point in the :py:attr:`brain.raw_data` array

	:return: :py:attr:`brain.df`
	:rtype: Dataframe with four columns: x,y,z,value

.. py:method:: brain.plot_projections(df)

	Plots the x, y, and z projections of the input dataframe in a matplotlib plot

	:param pd.DataFrame df: Dataframe with columns: 'x','y','z'
	:returns: Matplotlib figure with three labeled scatterplots

.. py:method:: brain.find_distance(t,point)

	Find euclidean distance between a point on the line defined by t and a data point

	:param float t: float value defining point on the line
	:param array point: array [x,y,z] defining data point
	:returns: distance between the two points
	:rtype: float

.. py:method:: brain.find_min_distance(row)

	Find the point on the curve that produces the minimum distance between the point and the data point using scipy.optimize.minimize(:py:func:`brain.find_distance`)

	:param pd.Series row: row from dataframe in the form of a pandas Series
	:returns: point in the curve (xc, yc, zc) and r
	:rtype: floats

.. py:method:: brain.find_alpha(xc,yc,zc)

	Find alpha which is the angle that specifies the position of the point along the curve

	:param float xc: x position of closest point on curve to datapoint
	:param float yc: y position of closest point on curve to datapoint
	:param float zc: zposition of closest point on curve to datapoint
	:returns: alpha, angle along the curve
	:rtype: float

.. py:method:: brain.integrand(x)

	Function to integrate to calculate arclength

	:param float x: integer value for x
	:returns: arclength value for integrating
	:rtype: float

.. py:method:: brain.find_length(xc)

	:param float row: Postion in the x axis along the curve
	:returns: Length of the arc along the curve between the row and the vertex
	:rtype: float

.. py:method:: brain.dist_to_plane(xz,row)

	Find shortest distance between point and the plane

	:param list xz: List of form [x position, y position]
	:param pd.Series row: row from dataframe in the form of a pandas Series
	:returns: Distance between the specified point and the plane
	:rtype: float

.. py:method:: brain.find_theta(row,r,zc)

	Calculate theta for a row containing data point in relationship to the flat plane

	:param pd.Series row: row from dataframe in the form of a pandas Series
	:param float r: Shortest distance between the point and the math model
	:param float zc: Z position of the closest point in the curve to the data point
	:returns: theta, angle between point and the flat plane
	:rtype: float

.. py:method:: brain.calc_coord(row)

	Calculate alpah, r, theta for a particular row

	:param pd.Series row: row from dataframe in the form of a pandas Series
	:returns: pd.Series populated with coordinate of closest point on the math model, r, theta, and ac (arclength)
	:rtype: pd.Series 

.. py:method:: transform_coordinates()

	Transform coordinate system so that each point is defined relative to math model by (alpha,theta,r) (only applied to df_thresh

	:returns: appends columns r, xc, yc, zc, ac, theta to :py:attr:`brain.df_thresh`

.. py:method:: brain.subset_data(sample_frac)

	Takes a random sample of the data based on the value between 0 and 1 defined for sample_frac

	:param sample_frac: Value between 0 and 1 specifying proportion of the dataset that should be randomly sampled for plotting
	:type: float or none
	:returns: :py:attr:`brain.subset`

.. py:method:: brain.add_thresh_df(df)

	Adds dataframe of thresholded and transformed data to :py:attr:`brain.df_thresh`

	:param pd.DataFrame df: dataframe of thesholded and transformed data
	:returns: :py:attr:`brain.df_thresh`


.. py:class:: math_model(model)

	Object to contain attributes associated with the math model of a sample

	:param array model: Array of coefficients calculated by np.polyfit

	.. py:attribute:: math_model.cf

		Array of coefficients for the math model

	.. py:attribute:: math_model.p

		Poly1d function for the math model to allow calculation and plotting of the model


.. py:function:: process_sample(filepath)

	Process single sample through :py:class:`brain` class and saves df to csv

	:param str filepath: Complete filepath to h5 data file
	:returns: Saves dataframe to csv with name of the original data file 