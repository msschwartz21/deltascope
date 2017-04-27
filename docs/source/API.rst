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

	.. py:attribute:: brain.subset

		Random sample of the dataframe from :py:attr:`brain.df_thresh`

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

.. py:method:: brain.find_theta(row,yc,zc)

	Calculate theta for a row containing data point in relationship to the flat plane

	:param pd.Series row: row from dataframe in the form of a pandas Series
	:param float r: Shortest distance between the point and the math model
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

.. py:method:: brain.plot_model(sample_frac=0.5,cmap='plt.cm.Greys')

	Plot two planes, line model, and percentage of points. Data is downsampled based on the value between 0 and 1 defined for sample_frac

	The returned plotly figure object can be most easily visualized using `plotly.offline.iplot(fig,filename='example')

	:param sample_frac: Value between 0 and 1 specifying proportion of the dataset that should be randomly sampled for plotting
	:type: float or none
	:returns: Plotly figure object



.. py:class:: plane(model,xx,yy,zz)

	Class to contain model and data associated with a plane

	:param model: OLS fitted model 
	:param array xx: Meshgrid array of x dimension
	:param array yy: Meshgrid array of y dimension
	:param array zz: Meshgrid array of z dimension

	.. py:attribute:: plane.model

		OLS fitted model

	.. py:attribute:: plane.xx

		Meshgrid array of x dimension

	.. py:attribute:: plane.yy

		Meshgrid array of y dimension

	.. py:attribute:: plane.zz

		Meshgrid array of z dimension


.. py:class:: math_model(coef)

	Class to contain attribues and data associated with math model

	:param dict coef: Dictionary containing coefficients to define equation of math model
	:param dict p: Dictionary containing calculated coefficients for y and z parabola
	:param array x: Array containing x coordinates
	:param array y: Array containing y coordinates
	:param array z: Array containing z coordinates

	.. py:attribute:: math_model.coef

		Dictionary containing coefficients of each term of math model such that:

		.. math::

			y = ex + fz + g
			z = ax^2 + bx + cy + d

	.. py:attribute:: math_model.p

		Dictionary containing coefficients of terms in math model such that:

		.. math::

			y = ay*x^2 + by*x + cy
			z = az*x^2 + bz*x + cz

	.. py:attribute:: math_model.x 

		Array containing x coordinates of model

	.. py:attribute:: math_model.y

		Array containing y coordinates of model

	.. py:attribute:: math_model.z

		Array containing z coordinates of model

	.. py:attribute:: math_model.vx

		x position of the vertex

	.. py:attribute:: math_model.vy

		y position of the vertex

	.. py:attribute:: math_model.vz

		z position of the vertex

	.. py:attribute:: math_model.fx

		x position of the focus

	.. py:attribute:: math_model.fy

		y position of the focus

	.. py:attribute:: math_model.fz

		z position of the focus

.. py:method:: math_model.calc_y(t)

	Calculate y value according to a given t

	:param float t: t value along the curve
	:returns: y value corresponding to t
	:rtype: float

.. py:method:: math_model.calc_z(t)

	Calculate z value according to a gien t

	:param float t: t value along the curve
	:returns: z value corresponding to t
	:rtype: float

.. py:method:: math_model.find_vertex()

	Calculates the position of the vertex

	:returns: :py:attr:`math_model.vx`, :py:attr:`math_model.vy`, :py:attr:`math_model.vz`

.. py:method:: math_model.find_focus()

	Calculates the position of the focus

	:returns: :py:attr:`math_model.fx`, :py:attr:`math_model.fy`, :py:attr:`math_model.fz`


.. py:function:: process_sample(filepath)

	Process single sample through :py:class:`brain` class and saves df to csv

	:param str filepath: Complete filepath to h5 data file
	:returns: Saves dataframe to csv with name of the original data file 