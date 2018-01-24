import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from sklearn.preprocessing import normalize.scale
from scipy.stats import sem

class landmarks:
	'''
	Class to handle calculation of landmarks to describe structural data

	:param list percbins: (or None) Must be a list of integers between 0 and 100
	:param int rnull: (or None) When the r value cannot be calculated it will be set to this value

	.. py:attribute:: brain.lm_wt_rf

		pd.DataFrame, which wildtype landmarks will be added to

	.. py:attribute:: brain.lm_mt_rf

		pd.DataFrame, which mutant landmarks will be added to

	.. py:attribute:: brain.rnull

		Integer specifying the value which null landmark calculations will be set to

	.. py:attribute:: brain.percbins

		Integer specifying the percentiles which will be used to calculate landmarks
	'''

	def __init__(self,percbins=[10,50,90],rnull=15):

		self.lm_wt_rf = pd.DataFrame()
		self.lm_mt_rf = pd.DataFrame()

		self.rnull = rnull
		self.percbins = percbins

	def calc_bins(self,Ldf,ac_num,tstep):
		'''
		Calculates alpha and theta bins based on ac_num and tstep

		Creates :py:attr:`landmarks.acbins` and :py:attr:`landmarks.tbins`

		.. warning:: `tstep` does not handle scenarios where 2pi is not evenly divisible by tstep

		:param dict Ldf: Dict dataframes that are being used for the analysis
		:param int ac_num: Integer indicating the number of divisions that should be made along alpha
		:param float tstep: The size of each bin used for alpha

		.. py:attribute:: landmarks.acbins

			List containing the boundaries of each bin along alpha based on `ac_num`

		.. py:attribute:: landmarks.tbins

			List containing the boundaries of each bin along theta based on `tstep`
		'''

		Ldf = Ldf.values()

		#Find the minimum and maximum values of alpha in the dataset
		acmin,acmax = 0,0
		for df in Ldf:
			if df.ac.min() < acmin:
				acmin = df.ac.min()
			if df.ac.max() > acmax:
				acmax = df.ac.max()

		#Correct min and max values so that arclength is equal on each side
		if abs(acmin) > acmax:
			self.acbins = np.linspace(acmin,abs(acmin),ac_num)
		else:
			self.acbins = np.linspace(-acmax,acmax,ac_num)

		#Calculate tbins divisions based on tstep
		self.tbins = np.arange(-np.pi,np.pi+tstep,tstep)

	def calc_perc(self,df,snum,dtype,out):
		'''
		Calculate landmarks for a dataframe based on the bins and percentiles that have been previously defined

		:param pd.DataFrame df: Dataframe containing columns x,y,z,alpha,r,theta
		:param str snum: String containing a sample identifier that can be converted to an integer
		:param str dtype: String describing the sample group to which the sample belongs, e.g. control or experimental
		:returns: pd.DataFrame with new landmarks appended
		'''

		D = {'stype':dtype}

		#Go through each arclength bins
		for a in range(len(self.acbins)):
			if a+1 < len(self.acbins):
				arange = [self.acbins[a],self.acbins[a+1]]

				#Go through each theta bin
				for t in range(len(self.tbins)):
					if t+1 < len(self.tbins):
						trange = [self.tbins[t],self.tbins[t+1]]

						#Go through each percentile
						for p in self.percbins:
							d = df[(df.ac > arange[0]) & (df.ac < arange[1])]
							d = d[(d.theta > trange[0]) & (d.theta < trange[1])]

							#Try to calculate percentile, but set to null if it fails
							try:
								r = np.percentile(d.r,p)
								pts = d[d.r < r].count()['i']
							except:
								r = self.rnull
								pts = 0

							#Create name of column based on bins
							L = []
							for s in [arange[0],arange[1],trange[0],trange[1]]:
								L.append(str(np.around(s,decimals=2)))
							name = '_'.join(L)

							D[name+'_'+str(p)+'_pts'] = pts
							D[name+'_'+str(p)+'_r'] = r

		out = out.append(pd.Series(D,name=int(snum)))
		return(out)

	def calc_wt_reformat(self,df,snum):
		'''
		.. warning:: Deprecated function, but includes code pertaining to calculating point based data
		'''

		D = {'stype':'wildtype'}

		for a in range(len(self.acbins)):
			if a+1 < len(self.acbins):
				arange = [self.acbins[a],self.acbins[a+1]]

				for t in range(len(self.tbins)):
					if t+1 < len(self.tbins):
						trange = [self.tbins[t],self.tbins[t+1]]

						for p in self.percbins:
							d = df[(df.ac > arange[0]) & (df.ac < arange[1])]
							d = d[(d.theta > trange[0]) & (d.theta < trange[1])]

							try:
								r = np.percentile(d.r,p)
								pts = d[d.r < r].count()['i']
							except:
								r = self.rnull
								pts = 0

							L = []
							for s in [arange[0],arange[1],trange[0],trange[1]]:
								L.append(str(np.around(s,decimals=2)))
							name = '_'.join(L)

							D[name+'_'+str(p)+'_pts-pts'] = pts
							D[name+'_'+str(p)+'_perc-pts'] = pts
							D[name+'_'+str(p)+'_pts-r'] = r
							D[name+'_'+str(p)+'_perc-r'] = r

		self.lm_wt_rf = self.lm_wt_rf.append(pd.Series(D,name=int(snum)))

	def calc_mt_landmarks(self,df,snum,wt):
		'''
		.. warning:: Deprecated function, but attempted to calculate mutant landmarks based on the number of points found in the wildtype standard
		'''

		D = {'stype':'mutant'}

		for c in wt.columns:
			if len(c.split('_')) == 6:
				amn,amx,tmn,tmx,p,dtype = c.split('_')
				p = int(p)

				# print(c)
				# print(df.count()['i'])
				d = df[(df.ac > float(amn))&(df.ac < float(amx))]
				# print(d.count()['i'])
				d = d[(d.theta > float(tmn))&(d.theta < float(tmx))]
				# print(d.count()['i'])

				if dtype.split('-')[0] == 'perc':
					try:
						perc_r = np.percentile(d.r,p)
						perc_pts = d[d.r < perc_r]
					except:
						perc_r = self.rnull
						perc_pts = 0

					if dtype.split('-')[1] == 'r':
						D[c] = perc_r
					else:
						D[c] = perc_pts

				else:
					if dtype.split('-')[1] == 'pts':
						#Calculate precent of points by dividing pts by total pts
						p = wt[c].mean()/d.count()['i']
						try:
							pt_r = np.percentile(d.r,p)
						except:
							pt_r = self.rnull
						D[c] = pt_r
					else:
						D[c] = np.nan

		self.lm_mt_rf = self.lm_mt_rf.append(pd.Series(D,name=int(snum)))

def reformat_to_cart(df):
	'''
	Take a dataframe in which columns contain the bin parameters and convert to a cartesian coordinate system

	:param pd.DataFrame df: Dataframe containing columns with string names that contain the bin parameter
	:returns: pd.DataFrame with each landmark as a row and columns: x,y,z,r,r_std,t,pts
	'''

	ndf = pd.DataFrame()
	for c in df.columns:
		if len(c.split('_')) == 6:
			amn,amx,tmn,tmx,p,dtype = c.split('_')
			x = np.mean([float(amn),float(amx)])
			t = np.mean([float(tmn),float(tmx)])

			if dtype == 'r':
				r = np.mean(df[c])
				r_std = sem(df[c])
				y = np.sin(t)*r
				z = np.cos(t)*r

				pts = np.mean(df['_'.join([amn,amx,tmn,tmx,p,'pts'])])

				D = pd.Series({'x':x,'y':y,'z':z,'r':r,'r_sem':r_std,'t':t,'pts':pts})

				ndf = ndf.append(pd.Series(D),ignore_index=True)

	return(ndf)

def convert_to_arr(xarr,tarr,DT,mdf,Ldf=[]):
	'''
	Convert a pandas dataframe containing landmarks as columns and samples as rows into a 3D numpy array

	The columns of `mdf` determine which landmarks will be saved into the array. Any additional dataframes that need to be converted can be included in Ldf

	:param np.array xarr: Array containing all unique x values of landmarks in the dataset
	:param np.array tarr: Array containing all unique t values of landmarks in the dataset
	:param str DT: Either ``r`` or ``pts`` indicating which data type should be saved to the array
	:param pd.DataFrame mdf: Main landmark dataframe containing landmarks as columns and samples as rows
	:param list Ldf: List of additional pd.DataFrames that should also be converted to arrays
	:returns: Array of the main dataframe and list of arrays converted from Ldf
	'''
	marr = np.zeros((len(xarr),len(tarr),len(mdf.index)))
	xarr = np.round(xarr,2)
	tarr = np.round(tarr,2)
	Larr = []
	for df in Ldf:
		Larr.append(np.zeros((len(xarr),len(tarr),len(df.index))))

	for c in mdf.columns:
		if len(c.split('_')) == 6:
			amn,amx,tmn,tmx,p,dtype = c.split('_')
			x = float(amn)#np.mean([float(amn),float(amx)])
			t = float(tmn) #np.mean([float(tmn),float(tmx)])

			if dtype == DT:
				marr[np.where(xarr==x)[0],np.where(tarr==t)[0]] = mdf[c]
				for arr,df in zip(Larr,Ldf):
					arr[np.where(xarr==x)[0],np.where(tarr==t)[0]] = df[c]

	return(marr,Larr)

def calc_variance(anum,dfs):
	'''
	Calculate the variance between samples according to bin position and variance between adjacent bins

	:param int anum: Number of bins which the arclength axis should be divided into
	:param dict dfs: Dictionary of dfs which are going to be processed
	:returns: Two arrays: svar (anum,tnum) and bvar (anum*tnum,snum)
	:rtype: np.array
	'''

	#Set up bins
	lm = landmarks(percbins=[50],rnull=15)
	lm.calc_bins(dfs.values(),anum,np.pi/4)

	 #Calculate landmarks
	outlm = pd.DataFrame()
	for k in dfs.keys():
		outlm = lm.calc_perc(dfs[k],k,'wt',outlm)

	#Convert to arr for variance calculation
	lmarr,arr = convert_to_arr(lm.acbins,lm.tbins,outlm)

	svar = np.var(lmarr,axis=2)

	Lvar = []
	for i in range(1,len(lm.acbins)-2):
		for t in range(0,len(lm.tbins)):
			#Save the variance of a particular bin and adjacent neighbors across a set of samples
			Lvar.append(np.var(lmarr[i-1:i+2,t],axis=0))

	return(svar,np.array(Lvar))

class anumSelect:

	def __init__(self,dfs):
		'''
		A class that assists in selecting the optimum value of anum

		:param dict dfs: Dictionary of pd.DataFrames with samples to use for optimization

		.. attribute:: anumSelect.dfs

			Dictionary of dataframes that will be used for the parameter sweep

		.. attribute:: anumSelect.Lsv

			List of sample variance arrays for each anum in the sweep

		.. attribute:: anumSelect.Lbv

			List of bin variance arrays for each anum in the sweep

		.. attribute:: anumSelect.Msv

			List of values of the average sample variance for each anum in the sweep

		.. attribute:: anumSelect.Mbv

			List of values of the average bin variance for each anum in the sweep

		.. attribute:: anumSelect.Llm

			List of landmark arrays for each anum in the sweep

		'''

		self.dfs = dfs
		self.Lsv,self.Lbv = [],[]
		self.Msv,self.Mbv = [],[]
		self.Llm = []

	def calc_variance(self,anum,tstep,percbins,rnull,DT):
		'''
		Calculate the variance between samples according to bin position and variance between adjacent bins

		:param int anum: Number of bins which the arclength axis should be divided into
		:param float tstep: The size of each bin used for alpha
		:param dict dfs: Dictionary of dfs which are going to be processed
		:param list percbins: Must be a list of integers between 0 and 100
		:param int rnull: When the r value cannot be calculated it will be set to this value
		:param str DT: Data type for which variance is measured, e.g. ``r`` or ``pts``
		'''

		#Set up bins
		lm = landmarks(percbins=percbins,rnull=rnull)
		lm.calc_bins(self.dfs,anum,tstep)

		#Calculate landmarks
		outlm = pd.DataFrame()
		for k in self.dfs.keys():
			outlm = lm.calc_perc(self.dfs[k],k,'s',outlm)

		#Convert to arr for variance calculations
		lmarr,arr = convert_to_arr(lm.acbins,lm.tbins,DT,outlm)
		self.Llm.append(lmarr)

		#Calculate variance between samples
		svar = np.nanvar(lmarr,axis=2)

		#Calculate variance between bins
		Lvar = []
		for i in range(1,len(lm.acbins-2)):
			for t in range(0,len(lm.tbins)):
				Lvar.append(np.nanvar(lmarr[i-1:i+2,t],axis=0))

		self.Lsv.append(svar)
		self.Msv.append(np.nanmean(svar))
		self.Lbv.append(np.array(Lvar))
		self.Mbv.append(np.nanmean(np.array(Lvar)))

		print(anum,'calculation complete')

	def param_sweep(self,tstep,amn=2,amx=50,astep=1,percbins=[50],rnull=15,DT='pts'):
		'''
		Calculate landmarks for each value of anum specified in input range

		:param float tstep: The size of each theta wedge in radians
		:param int amn: The minimum number of alpha bins that should be considered
		:param int amx: The maximum number of alpha bins that should be considered
		:param int astep: The step size in the range of amn to amx
		:param list percbins: (or None) Must be a list of integers between 0 and 100
		:param int rnull: (or None) When the r value cannot be calculated it will be set to this value
		:param str DT: Default=``pts`` Data type for which variance is measured, e.g. ``r`` or ``pts``

		.. attribute:: anumSelect.amn

			User defined minimum number of alpha bins considered in sweep

		.. attribute:: anumSelect.amx

			User defined maximum number of alpha bins considered in sweep

		.. attribute:: anumSelect.astep

			The step size in the range of amn to amx
		'''
		self.amn,self.amx,self.astep = amn,amx,astep

		for a in np.arange(self.amn,self.amx,self.astep):
			self.calc_variance(a,tstep,percbins,rnull,DT)

		print('Parameter sweep complete')

	def plot_rawdata(self):
		'''
		Plot raw data from parameter sweep
		'''
		x = np.arange(self.amn,self.amx,self.astep)
		fig,ax = plt.subplots()

		ax.plot(x,self.Mbv[2:],c='b',ls='--',label='Raw Bin Variance')
		ax.plot(x,self.Msv[2:],c='g',ls='--',label='Raw Sample Variance')

		ax.legend()
		ax.set_xlabel('Number of Alpha Bins')
		ax.set_ylabel('Relative Variance')

	def plot_fitted(self,dof):
		'''
		Plot polynomial fits over the raw data

		:param int dof: Degree of the polynomial function to be fit
		'''

		x = np.arange(self.amn,self.amx,self.astep)
		fig,ax = plt.subplots()

		ax.plot(x,self.Mbv[2:],c='b',ls='--',label='Raw Bin Variance')
		ax.plot(x,self.Msv[2:],c='g',ls='--',label='Raw Sample Variance')

		pbv = np.polyfit(x,self.Mbv[2:],dof)
		fbv = np.poly1d(pbv)
		ax.plot(x,fbv(x),c='b',label='Fit Bin Variance')

		psv = np.polyfit(x,self.Msv[2:],dof)
		fsv = np.poly1d(psv)
		ax.plot(x,fsv(x),c='g',label='Fit Sample Variance')

		ax.legend()
		ax.set_xlabel('Number of Alpha Bins')
		ax.set_ylabel('Relative Variance')

	def find_optimum_anum(self,dof,guess):
		'''
		Calculate optimum anum and plot

		:param int dof: Degree of the polynomial function to be fit
		:param int guess: Best guess of the optimum anum, which cannot be less than the maximum bin variance
		'''

		x = np.arange(self.amn,self.amx,self.astep)
		fig,ax = plt.subplots()

		pbv = np.polyfit(x,normalize(self.Mbv[2:])[0],dof)
		fbv = np.poly1d(pbv)
		ax.plot(x,fbv(x),c='b',label='Bin Variance')

		psv = np.polyfit(x,normalize(self.Msv[2:])[0],dof)
		fsv = np.poly1d(psv)
		ax.plot(x,fsv(x),c='g',label='Sample Variance')

		opt = minimize(fbv+fsv,guess)
		ax.axvline(opt.x,c='r',label='Optimum: '+str(np.round(opt.x[0],2)))

		ax.legend()
		ax.set_xlabel('Number of Alpha Bins')
		ax.set_ylabel('Relative Variance')

class graphSet:

	def __init__(self,tpairs,xarr,tarr):
		'''
		This object manages subobjects containing sample data, :class:`graphData`,
		and generates the appropriate landmark graph

		:param np.array tpairs: A list that specifies which theta bins should be paired together for graphing
		:param np.array xarr: An array listing the bin division points along alpha
		:param np.array tarr: Array listing the bin division points along theta

		.. attribute:: graphSet.tpairs

			A list that specifies which theta bins should be paired together for graphing

		.. attribute:: graphSet.xarr

			An array listing the bin division points along alpha

		.. attribute:: graphSet.tarr

			Array listing the bin division points along theta
		'''

		self.tpairs,self.xarr,self.tarr = tpairs,xarr,tarr

		self.Ls,self.Lc = [],[]
		self.Ds, self.Dc = {},{}

	def add_data(self,gD,stype,ctype,dtype):
		'''
		Add :class:`graphData` object with descriptive tags of stype,ctype,color

		:param obj gD: graphData object initialized for one sample set
		:param str stype: String describing the sample type
		:param str ctype: String describing channel type
		:param str dtype: Either ``pts`` or ``r`` specifies which type of landmark data will be used
		'''
		gD.prepare_data(self.xarr,self.tarr,dtype)

		self.add_to_dict(self.Ds,stype,ctype,gD)
		self.add_to_dict(self.Dc,ctype,stype,gD)

		#self.Ds.setdefault(stype).append({ctype:gD})
		#self.Dc.setdefault(ctype).append({stype:gD})
		self.Ls.append(stype)
		self.Lc.append(ctype)

	def add_to_dict(self,D,k1,k2,item):
		'''Add :class:`graphData` object to a dictionary checking for exisiting keys

		:param dict D: the dictionary which the data will be added to
		:param str k1: Key for the first index into the Dictionary
		:param str k2: Key for the second internal dictionary
		:param obj item: :class:`graphData` object
		'''

		if k1 in D.keys():
			D[k1][k2] = item
		else:
			D[k1] = {k2:item}

	def make_figure(self,a,figsize=(10,8),p=True):
		'''
		Creates a figure showing four theta slices and as many columns as ctypes

		.. todo:: P value scatter plot is broken

		.. todo:: Control p val for multiple testing

		:param float a: Alpha value for fill_between ribbons
		:param tuple figsize: Tuple specifying the height and width of the figure
		:param bool p: True if pvalue should be plotted

		.. attribute:: graphSet.fig

			Figure object created by :func:`graphSet.make_figure`

		.. attribute:: graphSet.axr

			Subplot axis array created by :func:`graphSet.make_figure`
		'''

		LsUn = np.unique(self.Ls)
		LcUn = np.unique(self.Lc)

		self.fig,self.axr = plt.subplots(4,len(LsUn),figsize=figsize,sharey=True)

		for j,c in enumerate(LcUn):
			dc = self.Dc[c]
			parr = stats.ttest_ind(dc[LsUn[0]].arr,dc[LsUn[1]].arr,axis=2,nan_policy='omit')[1]
			for i,p in enumerate(self.tpairs):
				for s in LsUn:
					go = dc[s]

					ti1 = np.where(self.tarr==p[0])[0][0]
					ti2 = np.where(self.tarr==p[1])[0][0]

					self.axr[i,j].fill_between(self.xarr,go.avg[:,ti1]+go.sem[:,ti1],go.avg[:,ti1]-go.sem[:,ti1],alpha=a,color=go.c,zorder=1)
					self.axr[i,j].fill_between(self.xarr,-go.avg[:,ti2]+go.sem[:,ti2],-go.avg[:,ti2]-go.sem[:,ti2],alpha=a,color=go.c,zorder=1)

					self.axr[i,j].plot(self.xarr,go.avg[:,ti1],c=go.c,zorder=2,label=c+s)
					self.axr[i,j].plot(self.xarr,-go.avg[:,ti2],c=go.c,zorder=2)

					if (s == 'mt') & (P==True):
						self.axr[i,j].scatter(self.xarr,go.avg[:,ti1],c=parr[:,ti1],cmap='Greys_r',zorder=3)
						self.axr[i,j].scatter(self.xarr,-go.avg[:,ti2],c=parr[:,ti2],cmap='Greys_r',zorder=3)
						print('plot pval')

				self.axr[i,j].legend()

class graphData:

	def __init__(self,rawdf,color):
		'''
		Object to contain data and attributes needed for graphing landmark dataset

		:param pd.DataFrame raw: Landmark dataframe with each column as a landmarks
		:param str color: Hexcode or single letter color code

		.. attribute:: graphData.c

			Color code for this dataset

		.. attribute:: graphData.rawdf

			Dataframe containing each landmark as a column
		'''
		self.c = color
		self.rawdf = rawdf

	def prepare_data(self,xarr,tarr,dtype):
		'''
		Convert to array and calculate average and sem

		:param arr xarr: List of min and max borders of the alpha bins
		:param arr tarr: List of min and max borders of the theta bins
		:param str dtype: Either ``pts`` or ``r`` specifies which data to save to array

		.. attribute:: graphData.arr

			Array containing the landmark data for this particular sample

		.. attribute:: graphData.avg

			Average of :attr:`graphData.arr`

		.. attribute:: graphData.sem

			Standard error of the mean for :attr:`graphData.arr`
		'''
		self.arr,L = convert_to_arr(xarr,tarr,dtype,self.rawdf)
		self.avg = np.nanmean(self.arr,axis=2)
		self.sem = stats.sem(self.arr,axis=2,nan_policy='omit')
