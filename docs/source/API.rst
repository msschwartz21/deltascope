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

	.. py:attribute:: brain.df_scl

		Dataframe containing data from :py:attr:`brain.df_thresh` after a scaling value has been applied

	.. py:attribute:: brain.scale

		Array with three values representing the constant by which to multiply x,y,z respectively

	.. py:attribute:: brain.pca

		PCA object managing the transformation matrix and any resulting transformations

	.. py:attribute:: brain.flip

		Boolean value indicating whether or not the sample needs to be inverted along the y axis

	.. py:attribute:: brain.vertex

		Array of shape [vx,vy,vz] containing values to use for translating sample

	.. py:attribute:: brain.df_align

		Dataframe containing point data aligned using PCA

	.. py:attribute:: brain.mm

		Math model object fit to data in brain object

.. py:method:: brain.read_data(filepath)

	Reads 3D data from file and selects appropriate channel based on the assumption that the channel with the most zeros has zero as the value for no signal

	:param str filepath: Filepath to hdf5 probability file
	:return: Array of shape [z,y,x] containing raw probability data

.. py:method:: brain.create_dataframe(data)

	Creates a pandas dataframe containing the x,y,z and signal/probability value for each point in the :py:attr:`brain.raw_data` array

	:param array data: Raw probability data in 3D array
	:return: Pandas DataFrame with xyz and probability value for each point

.. py:method:: brain.plot_projections(df,subset)

	Plots the x, y, and z projections of the input dataframe in a matplotlib plot

	:param pd.DataFrame df: Dataframe with columns: 'x','y','z'
	:param float subset: Value between 0 and 1 indicating what percentage of the df to subsample
	:returns: Matplotlib figure with three labeled scatterplots

.. py:method:: brain.preprocess_data(threshold,scale)

	Thresholds and scales data prior to PCA

	Creates :py:attr:`brain.threshold`, :py:attr:`brain.df_thresh`, and :py:attr:`brain.df_scl`

	:param float threshold: Value between 0 and 1 to use as a cutoff for minimum pixel value
	:param array scale: Array with three values representing the constant by which to multiply x,y,z respectively

.. py:method:: brain.calculate_pca()

	Create pca object and calculate transformation matrix in :py:attr:`brain.pca`

.. py:method:: brain.add_pca(pca)

	Add pca object from another channel and save as :py:attr:`brain.pca`

	:param sklearn.decomposition.PCA pca: PCA object containing transformation components that are already calculated

.. py:method:: brain.process_alignment_data(data,threshold,radius)

	Applies a double median filter to data to use for alignment

	:param array data: 3D array containing raw probability data
	:param float threshold: Value indicating the lower cutoff for positive signal
	:param int radius: Radius of the neighborhood that should be considered for the median filter
	:returns: dataframe after applying median filter twice and thresholding once
	:rtype: Pandas DataFrame

.. py:method:: brain.calculate_pca_median(data,threshold,radius)

	Calculate PCA transformation matrix based on data after applying median filter and threshold

	:param array data: 3D array containing raw probability data
	:param float threshold: Value between 0 and 1 indicating the lower cutoff for positive signal
	:param int radius: Radius of neighborhood that should be considered for the median filter

.. py:method:: brain.align_data(df,pca,comp_order,fit_dim[,deg=2,mm=None,vertex=None])

	Apply PCA transformation matrix and align data so that the vertex is at the origin

	:param pd.DataFrame df: dataframe containing thresholded xyz data
	:param array comp_order: Array specifies the assignment of components to x,y,z. Form [x component index, y component index, z component index], e.g. [0,2,1]
	:param array fit_dim: Array of length two containing two strings describing the first and second axis for fitting the model, e.g. ['x','z']
	:param deg: Degree of the function that should be fit to the model. deg=2 by default
	:type: int or None
	:param mm: Math model for primary channel
	:type: :py:class:`math_model` (:py:attr:`brain.mm`) or None
	:param vertex: Array indicating the translation values
	:type: Array [vx,vy,vz] (:py:attr:`brain.vertex`) or None

.. py:method:: brain.pca_transform(comp_order,fit_dim,flip_dim,[deg=2,mm=None,flip=None,vertex=None])

	Transforms data according to PCA transformation matrix and translates the sample so that the vertex of the parabola is at the origin

	Creates :py:attr:`brain.df_align`

	:param array comp_order: Array specifies the assignment of components to x,y,z. Form [x component index, y component index, z component index], e.g. [0,2,1]
	:param array fit_dim: Array of length two containing two strings describing the first and second axis for fitting the model, e.g. ['x','z']
	:param str flip_dim: String specifying the dimension in which to flip the data if necessary, e.g. 'z'
	:param deg: Degree of the function that should be fit to the model. deg=2 by default
	:type: int or None
	:param mm: Math model for primary channel
	:type: :py:class:`math_model` (:py:attr:`brain.mm`) or None
	:param flip: Boolean value to indicate whether or not the points need to be inverted along the y axis
	:type: Boolean (:py:attr:`brain.flip`) or None
	:param vertex: Array indicating the translation values
	:type: Array [vx,vy,vz] (:py:attr:`brain.vertex`) or None

.. py:method:: brain.fit_model(df,deg,fit_dim)

	Fit model to dataframe

	:param pd.DataFrame df: Dataframe containing at least x,y,z
	:param int deg: Degree of the function that should be fit to the model
	:param array fit_dim: Array of length two containing two strings describing the first and second axis for fitting the model, e.g. ['x','z']
	:returns: math model
	:rtype: :py:class:`math_model`

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

.. py:method:: embryo.process_channels(mthresh,gthresh,radius,scale,deg,primary_key,comp_order,fit_dim)
	
	Process all channels through the production of the :py:attr:`brain.df_align` dataframe

	:param float mthresh: Value between 0 and 1 to use as a cutoff for minimum pixel value for median data
	:param float gthresh: Value between 0 and 1 to use as a cutoff for minimum pixel value for general data
	:param int radius: Size of the neighborhood area to examine with median filter
	:param array scale: Array with three values representing the constant by which to multiply x,y,z respectively
	:param int deg: Degree of the function that should be fit to the model
	:param str primary_key: Key for the primary structural channel which PCA and the model should be fit too
	:param array comp_order: Array specifies the assignment of components to x,y,z. Form [x component index, y component index, z component index], e.g. [0,2,1]
	:param array fit_dim: Array of length two containing two strings describing the first and second axis for fitting the model, e.g. ['x','z']

.. py:method:: embryo.save_projections(subset)

	Save projections of both channels into png files in :py:attr:`embryo.outdir` following the naming scheme [:py:attr:`embryo.name`]_[:py:attr:`embryo.number`]_[`channel name`]_MIP.png

	:param float subset: Value between 0 and 1 to specify the fraction of the data to randomly sample for plotting

.. py:method:: embryo.save_psi()

	Save all channels into psi files following the naming scheme [:py:attr:`embryo.name`]_[:py:attr:`embryo.number`]_[`channel name`].psi

.. py:method:: embryo.add_psi_data(filepath,key)

	Read psi data into a channel dataframe

	:param str filepath: Complete filepath to data
	:param str key: Descriptive key for channel dataframe in dictionary


.. py:class:: math_model(model)

	Object to contain attributes associated with the math model of a sample

	:param array model: Array of coefficients calculated by np.polyfit

	.. py:attribute:: math_model.cf

		Array of coefficients for the math model

	.. py:attribute:: math_model.p

		Poly1d function for the math model to allow calculation and plotting of the model


.. py:function:: process_sample(num,root,outdir,name,chs,prefixes,threshold,scale,deg,primary_key)

	Process single sample through :py:class:`brain` class and saves df to csv

	:param str num: Sample number
	:param str root: Complete path to the root directory for this sample set
	:param str name: Name describing this sample set
	:param str outdir: Complete path to output directory
	:param array chs: Array containing strings specifying the directories for each channel
	:param array prefixes: Array containing strings specifying the file prefix for each channel
	:param float threshold: Value between 0 and 1 to use as a cutoff for minimum pixel value
	:param array scale: Array with three values representing the constant by which to multiply x,y,z respectively
	:param int deg: Degree of the function that should be fit to the model
	:param str primary_key: Key for the primary structural channel which PCA and the model should be fit too



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

.. py:function:: calculate_models(Ldf)

	Calculate model for each dataframe in list and add to new dataframe

	:param list Ldf: List of dataframes containing aligned data
	:returns: Dataframe with a,b,c values for parabolic model
	:rtype: pd.DataFrame

.. py:function:: fit_distribution(y,var,distname)

	Fit specified distribution to downsample data and return parameters and figure with plotted data

	:param pd.DataFrame y: Dataframe with complete data
	:param str var: Name of column in variable, which needs to be fitting
	:param str distname: 'beta' or 'uniform'
	:returns: Matplotlib figure object and parameter list

.. py:function:: calculate_rms(df)

	Calculate root mean squared error for a sample using r values

	:param pd.DataFrame df: Dataframe containing 'r'
	:returns: RMSE value
	:rtype: float

.. py:function:: test_distributions(y,snum[,plot=False,outdir=None])

	Tests fit of gamma, beta, and normal distributions on preprocessed data

	:param array y: Array of preprocessed data to be used for distribution fit
	:param str snum: String with sample number
	:param plot: Set to True to receive plotted data in outdir
	:type: Boolean or None
	:param str outdir: Set to path to outdir otherwise files will save in current working directory
	:returns: Dataframe containing parameters and KS test output for each distribution

.. py:function:: test_beta(y,snum[,plot=False,outdir=None])

	Fits beta distrubtion to data and returns parameters in df

	:param array y: Array of preprocessed data to be used for distribution fit
	:param str snum: String with sample number
	:param plot: Set to True to receive plotted data in outdir
	:type: Boolean or None
	:param str outdir: Set to path to outdir otherwise files will save in current working directory
	:returns: Dataframe with parameters and KS test result
	:rtype: pd.DataFrame

.. py:function:: test_gamma(y,snum[,plot=False,outdir=None])

	Fits gamma distrubtion to data and returns parameters in df

	:param array y: Array of preprocessed data to be used for distribution fit
	:param str snum: String with sample number
	:param plot: Set to True to receive plotted data in outdir
	:type: Boolean or None
	:param str outdir: Set to path to outdir otherwise files will save in current working directory
	:returns: Dataframe with parameters and KS test result
	:rtype: pd.DataFrame

.. py:function:: fit_gamma(y[,plot=False])

	Fit a gamma distribution to data and return parameter and plot if specified

	:param array y: Array of preprocessed data to be used for distribution fit
	:param plot: Set to True to receive plotted data in outdir
	:type: Boolean or None

.. py:function:: fit_beta(y[,plot=False])

	Fit a beta distribution to data and return parameter and plot if specified

	:param array y: Array of preprocessed data to be used for distribution fit
	:param plot: Set to True to receive plotted data in outdir
	:type: Boolean or None

.. py:function:: fit_norm(y[,plot=False])

	Fit a normal distribution to data and return parameter and plot if specified

	:param array y: Array of preprocessed data to be used for distribution fit
	:param plot: Set to True to receive plotted data in outdir
	:type: Boolean or None

.. py:function:: calculate_sample_error(y,param,distname,n)

	Calculate RMSE based on KDE of a sample and PDF of a distribution

	:param array y: Array of preprocessed data to be used for distribution fit
	:param array param: Array of parameters for specified distribution
	:param str distname: 'beta' or 'gamma'
	:param int n: Length of x for fitting data
	:returns: Root mean squared error value
	:rtype: float

.. py:function:: beta_pdf(x,param)

	Fit beta distribution to array x to get probability distribution function

	:param array x: Array to fit pdf to
	:param array param: Array of [a,b,loc,scale] for beta distribution
	:returns: Beta probability distribution function
	:rtype: array

.. py:function:: gamma_pdf(x,param)

	Fit gamma distribution to array x to get probability distribution function

	:param array x: Array to fit pdf to
	:param array param: Array of [a,loc,scale] for gamma distribution
	:returns: Gamma probability distribution function
	:rtype: array

.. py:function:: norm_pdf(x,param)

	Fit normal distribution to array x to get probability distribution function

	:param array x: Array to fit pdf to
	:param array param: Array of [loc,scale] for normal distribution
	:returns: Normal probability distribution function
	:rtype: array

.. py:function:: beta_rvs(param)

	Generate sampling of random points from specified beta distribution

	:param array param: Array of [a,b,loc,scale] for beta distribution
	:returns: Array of random values from distribution
	:rtype: array

.. py:function:: gamma_rvs(param)

	Generate sampling of random points from specified gamma distribution

	:param array param: Array of [a,loc,scale] for gamma distribution
	:returns: Array of random values from distribution
	:rtype: array

.. py:function:: calculate_xcurve

	.. warning:: not working

.. py:function:: calculate_zcurve

	.. warning:: not working

.. py:function:: calculate_y

	.. warning:: not working

.. py:function:: calculate_xyz

	.. warning:: not working

.. py:function:: generate_poc

	.. warning:: not working

.. py:function:: read_psi_to_dict(directory,dtype)

	Read psis from directory into dictionary of dfs with filtering based on dtype

	:param str directory: Directory to get psis from 
	:param str dtype: Usually 'AT' or 'ZRF1'
	:returns: Dictionary of pd.DataFrame
	:rtypes: dictionary

.. py:function:: concatenate_dfs(dfdict)

	Concatenated a dictionary of dfs into one df

	:param dict dfdict: Dictionary of pd.DataFrames
	:returns: Concatenated df
	:rtype: pd.DataFrame

.. py:function:: generate_kde(data,var,x[,absv=False])

	Generate list of KDEs from either dictionary or list of data

	:param data: pd.DataFrames to convert
	:type: dict or list
	:param str var: Name of column to select from df
	:param array x: Array of datapoints to evaluate KDE on
	:param absv: Set to True to use absolute value of selected data for KDE calculation
	:type: boolean or None
	:returns: List of KDE arrays

.. py:function:: fit_bimodal_theta(D,split,frac,x)

	Fit two distributions to bimodal theta data from dict D

	:param dict D: dictionary of pd.DataFrame
	:param split: Value to divide data for two distributions
	:type: int or float
	:param float frac: Fraction of the data to sample for distribution fitting
	:param array x: Array of datapoint to evaluate pdf on 
	:returns: Single pdf array containing both distributions scaled to proportion of points in each half of data
	:rtype: array

.. py:function:: calculate_area_error(pdf,Lkde,x)

	Calculate area between PDF and each kde in Lkde

	:param array pdf: Array of probability distribution function that is the same shape as kdes in Lkde
	:param list Lkde: List of arrays of Kdes 
	:param array x: Array of datapoints used to generate pdf and kdes
	:returns: List of error values for each kde in Lkde
	:rtype: list

.. py:function:: rescale_variable(Ddfs,var,newvar)

	Rescale variable from -1 to 1 and save in newvar column on original dataframe

	:param dict Ddfs: Dictionary of pd.DataFrames
	:param str var: Name of column to select from dfs
	:param str newvar: Name to use for new data in appended column
	:returns: Dictionary of dataframes containing column of rescaled data

