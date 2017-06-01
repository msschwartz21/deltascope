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

	.. py:attribute:: brain.subset

		Random sample of the dataframe from :py:attr:`brain.df_thresh`

	.. py:attribute:: brain.df_align

		Dataframe containing point data aligned using PCA

	.. py:attribute:: brain.mm

		Math model object fit to data in brain object

.. py:method:: brain.read_data(filepath)

	Reads 3D data from file and selects appropriate channel based on the assumption that the channel with the most zeros has zero as the value for no signal

	:param str filepath: Filepath to hdf5 probability file
	:return: Array of shape [z,y,x] containing raw probability data

.. py:method:: brain.create_dataframe()

	Creates a pandas dataframe containing the x,y,z and signal/probability value for each point in the :py:attr:`brain.raw_data` array

	:return: :py:attr:`brain.df`
	:rtype: Dataframe with four columns: x,y,z,value

.. py:method:: brain.plot_projections(df,subset)

	Plots the x, y, and z projections of the input dataframe in a matplotlib plot

	:param pd.DataFrame df: Dataframe with columns: 'x','y','z'
	:param float subset: Value between 0 and 1 indicating what percentage of the df to subsample
	:returns: Matplotlib figure with three labeled scatterplots

.. py:method:: brain.align_sample(threshold,scale[,deg=2])

	Realigns sample axes using PCA and translates so that the vertex is at the origin

	:param float threshold: Value between 0 and 1 to use as a cutoff for minimum pixel value
	:param array scale: Array with three values representing the constant by which to multiply x,y,z respectively
	:param deg: Degree of the function that should be fit to the model. deg=2 by default
	:type: int or None
	:returns: :py:attr:`brain.df_align` and :py:attr:`brain.mm`

.. py:method:: brain.fit_model(df,deg)

	Fit model to dataframe

	:param pd.DataFrame df: Dataframe containing at least x,y,z
	:param int deg: Degree of the function that should be fit to the model

.. py:method:: brain.find_distance(t,point)

	Find euclidean distance between math model(t) and data point in the xy plane

	:param float t: float value defining point on the line
	:param array point: array [x,y] defining data point
	:returns: distance between the two points
	:rtype: float

.. py:method:: brain.find_min_distance(row)

	Find the point on the curve that produces the minimum distance between the point and the data point using scipy.optimize.minimize(:py:func:`brain.find_distance`)

	:param pd.Series row: row from dataframe in the form of a pandas Series
	:returns: point in the curve (xc, yc, zc) and r
	:rtype: floats

.. py:method:: brain.integrand(x)

	Function to integrate to calculate arclength

	:param float x: integer value for x
	:returns: arclength value for integrating
	:rtype: float

.. py:method:: brain.find_arclength(xc)

	Calculate arclength by integrating the derivative of the math model in xy plane

	.. math:: 

		\int_{vertex}^{point} \sqrt{1 + (2ax + b)^2}

	:param float row: Postion in the x axis along the curve
	:returns: Length of the arc along the curve between the row and the vertex
	:rtype: float

.. py:method:: brain.find_theta(row,xc,zc)

	Calculate theta for a row containing data point in relationship to the xy plane

	:param pd.Series row: row from dataframe in the form of a pandas Series
	:param float xc: X position of the closest point in the curve to the data point
	:param float zc: Z position of the closest point in the curve to the data point
	:returns: theta, angle between point and the xy plane
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


.. py:class:: embryo(name,number,outdir)

	Class to managed multiple brain objects in a multichannel sample

	:param str name: Name of this sample set
	:param str number: Sample number corresponding to this embryo
	:param str outdir: Path to directory for output files

	.. py:attribute:: embryo.chnls

		Dictionary containing the :py:class:`brain` object for each channel

	.. py:attribute:: embryo.outdir

		Path to directory for output files

	.. py:attribute:: embryo.name

		Name of this sample set

	.. py:attribute:: embryo.number

		Sample number corresponding to this embryo

.. py:method:: embryo.add_channel(filepath,key)

	Add channel to :py:attr:`embryo.chnls` dictionary

	:param str filepath: Complete filepath to image
	:param str key: Name of the channel

.. py:method:: embryo.process_channels(threshold,scale,deg)
	
	Process all channels through the production of the :py:attr:`brain.df_align` dataframe

	:param float threshold: Value between 0 and 1 to use as a cutoff for minimum pixel value
	:param array scale: Array with three values representing the constant by which to multiply x,y,z respectively
	:param int deg: Degree of the function that should be fit to the model

.. py:method:: embryo.save_projections(subset)

	Save projections of both channels into png files in :py:attr:`embryo.outdir` following the naming scheme [:py:attr:`embryo.name`]_[:py:attr:`embryo.number`]_[`channel name`]_MIP.png

	:param float subset: Value between 0 and 1 to specify the fraction of the data to randomly sample for plotting

.. py:method:: embryo.save_psi()

	Save all channels into psi files following the naming scheme [:py:attr:`embryo.name`]_[:py:attr:`embryo.number`]_[`channel name`].psi



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



.. py:function:: write_header(f)

	Writes header for PSI file with columns Id,x,y,z,ac,r,theta

	:param file f: file object created by 'open(filename,'w')`

.. py:function:: write_data(filepath,df)

	Writes data in PSI format to file after writing header using :py:func:`write_header`. Closes file at the conclusion of writing data.

	:param str filepath: Complete filepath to output file
	:param pd.DataFrame df: dataframe containing columns x,y,z,ac,r,theta

.. py:function:: read_psi(filepath)

	Reads psi file at the given filepath and returns data in a pandas DataFrame

	:param str filepath: Complete filepath to file
	:returns: Dataframe containing data
	:rtype: pd.DataFrame