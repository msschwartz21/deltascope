import os
import numpy as np
import h5py
import skimage.io as io
import matplotlib.pyplot as plt
import time
import itertools as it
import matplotlib as mpl
import multiprocessing as mp
from functools import partial
import statsmodels.formula.api as smf
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
#import plotly.plotly as py
#import plotly.graph_objs as go
from scipy.optimize import minimize
import scipy
from sklearn.decomposition import PCA
from skimage.filters import median
from skimage.morphology import disk
from sklearn.metrics import mean_squared_error
from scipy.integrate import simps

class brain:

	def __init__(self):
		'''Initialize brain object'''

	def read_data(self,filepath):
		'''Reads 3D data from file and selects appropriate channel'''

		#Read h5 file and extract probability data
		f = h5py.File(filepath,'r')

		#Select either channe0/1 or exported_data
		d = f.get('exported_data')
		if d != None:
			c1 = np.array(d[:,:,:,0])
			c2 = np.array(d[:,:,:,1])
		elif d == None:
			c1 = np.array(f.get('channel0'))
			c2 = np.array(f.get('channel1'))

		#Figure out which channels has more zeros and therefore is background
		if np.count_nonzero(c1<0.1) > np.count_nonzero(c1>0.9):
			self.raw_data = c1
		else:
			self.raw_data = c2

	def create_dataframe(self,data,scale):
		'''Creates a pandas dataframe containing the x,y,z and signal/probability value 
		for each point in the :py:attr:`brain.raw_data` array'''

		dim = data.shape
		xyz = np.zeros((dim[0],dim[1],dim[2],4))

		#Generate array with xyz values for each point
		for x in range(dim[2]):

			xyz[:,:,x,2] = x
			xyz[:,:,x,3] = data[:,:,x]

			zindex = np.arange(0,dim[0],1)
			yindex = np.arange(0,dim[1],1)
			gy,gz = np.meshgrid(yindex,zindex)

			xyz[:,:,x,0] = gz
			xyz[:,:,x,1] = gy

		flat = np.reshape(xyz,(-1,4))

		#Create dataframe of points and scale to microns
		df = pd.DataFrame({'x':flat[:,2]*scale[0],'y':flat[:,1]*scale[1],'z':flat[:,0]*scale[2],'value':flat[:,3]})
		return(df)

	def plot_projections(self,df,subset):
		'''Plots x, y, and z projections from the dataframe'''

		df = df.sample(frac=subset)
    
    	#Create figure and subplots
		fig = plt.figure(figsize=(12,6))
		ax = fig.add_subplot(131)
		ay = fig.add_subplot(132)
		az = fig.add_subplot(133)

		#Create scatter plot for each projection
		ax.scatter(df.x,df.z)
		ay.scatter(df.x,df.y)
		az.scatter(df.z,df.y)

		#Plot model
		xvalues = np.arange(np.min(df.x),np.max(df.x))
		ax.plot(xvalues,self.mm.p(xvalues),c='y')

		#Add labels
		ax.set_title('Y projection')
		ay.set_title('Z projection')
		az.set_title('X projection')
		ax.set_xlabel('X')
		ax.set_ylabel('Z')
		ay.set_xlabel('X')
		ay.set_ylabel('Y')
		az.set_xlabel('Z')
		az.set_ylabel('Y')

		#Adjust spacing and show plot
		plt.subplots_adjust(wspace=0.4)
		
		return(fig)

	def preprocess_data(self,threshold,scale,microns):
		'''Thresholds and scales data prior to PCA'''

		self.df = self.create_dataframe(self.raw_data,microns)

		#Create new dataframe with values above threshold
		self.threshold = threshold
		self.df_thresh = self.df[self.df.value > self.threshold]

		#Scale xyz by value in scale array to force PCA axis selection
		self.scale = scale
		self.df_scl = pd.DataFrame({
			'x':self.df_thresh.x * self.scale[0], 
			'y':self.df_thresh.y * self.scale[1],
			'z':self.df_thresh.z * self.scale[2]})

	def process_alignment_data(self,data,threshold,radius,microns):
		'''Applies a double median filter to data to use for alignment'''

		#Iterate over each plane and apply median filter twice
		out = np.zeros(data.shape)
		for z in range(data.shape[0]):
			out[z] = median(median(data[z],disk(radius)),disk(radius))

		outdf = self.create_dataframe(out,microns)
		thresh = outdf[outdf.value > threshold]
		return(thresh)

	def calculate_pca_median(self,data,threshold,radius,microns):
		'''Transform data according to median filtered raw data'''

		self.median = self.process_alignment_data(data,threshold,radius,microns)

		self.pcamed = PCA()
		self.pcamed.fit(self.median[['x','y','z']])
	
	def align_data(self,df,pca,comp_order,fit_dim,deg=2,mm=None,vertex=None,flip=None):
		'''Apply PCA transformation matrix and align data so that 
		the vertex is at the origin'''

		#PCA transformation matrix
		fit = pca.transform(df[['x','y','z']])
		df_fit = pd.DataFrame({
			'x':fit[:,comp_order[0]],
			'y':fit[:,comp_order[1]],
			'z':fit[:,comp_order[2]]
			})

		#If vertex for translation is not included
		if vertex == None:
			#Calculate model
			model = np.polyfit(df_fit[fit_dim[0]],df_fit[fit_dim[1]],deg=deg)
			p = np.poly1d(model)

			#Find vertex
			a = -model[1]/(2*model[0])
			if fit_dim[0] == 'x':
				vx = a
				if fit_dim[1] == 'y':
					vy = p(vx)
					vz = df_fit.z.mean()
				else:
					vz = p(vx)
					vy = df_fit.y.mean()
			elif fit_dim[0] == 'y':
				vy = a
				if fit_dim[1] == 'x':
					vx = p(a)
					vz = df_fit.z.mean()
				else:
					vz = p(a)
					vx = df_fit.x.mean()
			elif fit_dim[0] == 'z':
				vz = a
				if fit_dim[1] == 'x':
					vx = p(a)
					vy = df_fit.y.mean()
				else:
					vy = p(a)
					vx = df_fit.x.mean()
			self.vertex = [vx,vy,vz]
		else:
			self.vertex = vertex

		#Translate data so that the vertex is at the origin
		self.df_align = pd.DataFrame({
			'x': df_fit.x - self.vertex[0],
			'y': df_fit.y - self.vertex[1],
			'z': df_fit.z - self.vertex[2]
			})

		#Rotate data by 180 degrees if necessary
		if flip == None:
			#Calculate model
			model = np.polyfit(df_fit[fit_dim[0]],df_fit[fit_dim[1]],deg=deg)
			p = np.poly1d(model)

			#If a is less than 0, rotate data
			if model[0] < 0:
				self.df_align = self.flip_data(self.df_align)
		elif flip == True:
			self.df_align = self.flip_data(self.df_align)

		#Calculate final math model
		if mm == None:
			self.mm = self.fit_model(self.df_align,deg,fit_dim)
		else:
			self.mm = mm

	def flip_data(self,df):
		'''Rotate data by 180 degrees'''

		r = np.array([[np.cos(np.pi),0,np.sin(np.pi)],
			[0,1,0],
			[-np.sin(np.pi),0,np.cos(np.pi)]])

		rot = np.dot(np.array(df),r)

		dfr = pd.DataFrame({'x':rot[:,0],'y':rot[:,1],'z':rot[:,2]})
		return(dfr)

	def fit_model(self,df,deg,fit_dim):
		'''Fit model to dataframe'''

		mm = math_model(np.polyfit(df[fit_dim[0]], df[fit_dim[1]], deg=deg))
		return(mm)

	##### Deprecated PCA functions

	def pca_transform(self,comp_order,fit_dim,flip_dim,deg=2,mm=None,flip=None,vertex=None):
		'''Transform data according to PCA fit parameters'''

		pca_fit = self.pca.transform(self.df_scl[['x','y','z']])

		#Create pca dataframe and remove scaling
		df_unscl = pd.DataFrame({
			'x':pca_fit[:,comp_order[0]]/self.scale[0],
			'y':pca_fit[:,comp_order[1]]/self.scale[1],
			'z':pca_fit[:,comp_order[2]]/self.scale[2]
			})

		#if flip = None, then this is the primary channel
		if flip == None:
			#Find first model
			model = np.polyfit(df_unscl[fit_dim[0]],df_unscl[fit_dim[1]],deg=deg)
			p = np.poly1d(model)

			#If parabola is upside down, flip y coordinates and set flip value
			if model[0] < 0:
				self.flip = True
				df_unscl[flip_dim] = df_unscl[flip_dim] * -1
				#Recalculate model
				model = np.polyfit(df_unscl[fit_dim[0]],df_unscl[fit_dim[1]],deg=deg)
				p = np.poly1d(model)
			else:
				self.flip = False
		else:
			if flip == True:
				df_unscl[flip_dim] = df_unscl[flip_dim] * -1

		#If vertex for translation is not included
		if vertex == None:
			#Find vertex
			a = -model[1]/(2*model[0])
			if fit_dim[0] == 'x':
				vx = a
				if fit_dim[1] == 'y':
					vy = p(vx)
					vz = df_unscl.y.mean()
				else:
					vz = p(vx)
					vy = df_unscl.z.mean()
			elif fit_dim[0] == 'y':
				vy = a
				if fit_dim[1] == 'x':
					vx = p(a)
					vz = df_unscl.z.mean()
				else:
					vz = p(a)
					vx = df_unscl.x.mean()
			elif fit_dim[0] == 'z':
				vz = a
				if fit_dim[1] == 'x':
					vx = p(a)
					vy = df_unscl.y.mean()
				else:
					vy = p(a)
					vx = df_unscl.x.mean()
			self.vertex = [vx,vy,vz]
		else:
			vx = vertex[0]
			vy = vertex[1]
			vz = vertex[2]

		#Translate data so that the vertex is at the origin
		self.df_align = pd.DataFrame({
			'x': df_unscl.x - vx,
			'y': df_unscl.y - vy,
			'z': df_unscl.z - vz
			})

		if mm == None:
			#Calculate final model based on translated and aligned data
			self.mm = self.fit_model(self.df_align,deg,fit_dim)
		else:
			self.mm = mm

	def calculate_pca(self):
		'''Create pca object and calculate transformation matrix'''

		self.pca = PCA(n_components=3)
		self.pca.fit(self.df_scl[['x','y','z']])

	def add_pca(self,pca):
		'''Add pca object from another channel'''

		self.pca = pca

	def pca_double_transform(self,threshold,fit_dim,pca1=None,pca2=None,deg=2,mm=None,vertex=None):
		'''Transform data twice to align in plane'''

		#Create new dataframe with values above threshold
		self.threshold = threshold
		self.df_thresh = self.df[self.df.value > self.threshold]

		#Calulate first transformation matrix
		if pca1 == None:
			self.pca1 = PCA(n_components=3)
			self.pca1.fit(self.df_thresh[['x','y','z']])
		else:
			self.pca1 = pca1

		#Transform data according to assigned components
		fit1 = self.pca1.transform(self.df_thresh[['x','y','z']])
		#Create array for second pca transformation
		a = np.zeros((fit1.shape[0],2))
		a[:,0] = fit1[:,1]
		a[:,1] = fit1[:,2]

		if pca2 == None:
			self.pca2 = PCA(n_components=2)
			self.pca2.fit(a)
		else:
			self.pca2 = pca2

		#Transform second according to second PCA
		fit2 = self.pca2.transform(a)

		#Put data in inbetween dataframe
		df = pd.DataFrame({
			'x': fit1[:,0],
			'y': fit2[:,1],
			'z': fit2[:,0]
			})

		if mm == None:
			model = self.fit_model(df,deg,fit_dim)
		else:
			self.mm = mm

		#If vertex for translation is not included
		if vertex == None:
			#Find vertex
			a = -model.cf[1]/(2*model.cf[0])
			if fit_dim[0] == 'x':
				vx = a
				if fit_dim[1] == 'y':
					vy = model.p(vx)
					vz = df.y.mean()
				else:
					vz = model.p(vx)
					vy = df.z.mean()
			elif fit_dim[0] == 'y':
				vy = a
				if fit_dim[1] == 'x':
					vx = model.p(a)
					vz = df.z.mean()
				else:
					vz = model.p(a)
					vx = df.x.mean()
			elif fit_dim[0] == 'z':
				vz = a
				if fit_dim[1] == 'x':
					vx = model.p(a)
					vy = df.y.mean()
				else:
					vy = model.p(a)
					vx = df.x.mean()
			self.vertex = [vx,vy,vz]
		else:
			vx = vertex[0]
			vy = vertex[1]
			vz = vertex[2]

		#Translate data so that the vertex is at the origin
		self.df_align = pd.DataFrame({
			'x': df.x - vx,
			'y': df.y - vy,
			'z': df.z - vz
			})

		if mm == None:
			#Calculate final model based on translated and aligned data
			self.mm = self.fit_model(self.df,deg,fit_dim)
		else:
			self.mm = mm
	###### Functions associated with alpha, r, theta coordinate system ######

	def find_distance(self,t,point):
		'''Find euclidean distance between math model(t) and data point in the xz plane'''

		x = float(t)
		z = self.mm.p(x)

		#Calculate distance between two points passed as array
		dist = np.linalg.norm(point - np.array([x,z]))

		return(dist)

	def find_min_distance(self,row):
		'''Find the point on the curve that produces the minimum distance between the point and the data point'''

		dpoint = np.array([row.x,row.z])

		#Use scipy.optimize.minimize to find minimum solution of brain.find_distance, 
		#dpoint[0] is starting guess for x value
		result = minimize(self.find_distance, dpoint[0], args=(dpoint))

		x = result['x'][0]
		y = 0
		z = self.mm.p(x)
		r = result['fun']

		return(x,y,z,r)
		#return(pd.Series({'xc':x, 'yc':y, 'zc':z, 'r':r}))

	def integrand(self,x):
		'''Function to integrate to calculate arclength between vertex and x'''

		y_prime = self.mm.cf[0]*2*x + self.mm.cf[1]

		arclength = np.sqrt(1 + y_prime**2)
		return(arclength)

	def find_arclength(self,xc):
		'''Calculate arclength for a row'''

		ac,err = scipy.integrate.quad(self.integrand,xc,0)
		return(ac)

	def find_theta(self,row,zc,yc):
		'''Find theta value for a row describing angle between point and plane'''

		theta = np.arctan2(row.y-yc,row.z-zc)
		return(theta)

	def calc_coord(self,row):
		'''Calculate alpah, r, theta for a particular row'''

		xc,yc,zc,r = self.find_min_distance(row)
		ac = self.find_arclength(xc)
		theta = self.find_theta(row,zc,yc)

		return(pd.Series({'xc':xc, 'yc':yc, 'zc':zc,
					'r':r, 'ac':ac, 'theta':theta}))

	def transform_coordinates(self):
		'''Transform coordinate system so that each point is defined relative to 
		math model by (alpha,theta,r) (only applied to df_align'''

		#Calculate alpha, theta, r for each row in dataset
		self.df_align = self.df_align.join(self.df_align.apply((lambda row: self.calc_coord(row)), axis=1))

	def subset_data(self,sample_frac=0.5):
		'''Subset data based on proportion set in sample_frac'''

		self.subset = self.df_align.sample(frac=sample_frac)

	def add_thresh_df(self,df):
		'''Add dataframe of thresholded and transformed data to self.df_thresh'''

		self.df_thresh = df

	def add_aligned_df(self,df):
		'''Add dataframe of aligned data'''

		self.df_align = df

class embryo:
	'''Class to managed multiple brain objects in a multichannel sample'''

	def __init__(self,name,number,outdir):
		'''Initialize embryo object'''

		self.chnls = {}
		self.outdir = outdir
		self.name = name
		self.number = number

	def add_channel(self,filepath,key):
		'''Add channel to self.chnls dictionary'''

		s = brain()
		s.read_data(filepath)

		self.chnls[key] = s

	def process_channels(self,mthresh,gthresh,radius,scale,microns,deg,primary_key,comp_order,fit_dim):
		'''Process channels through alignment'''

		#Process primary channel
		self.chnls[primary_key].preprocess_data(gthresh,scale,microns)

		self.chnls[primary_key].calculate_pca_median(self.chnls[primary_key].raw_data,
			mthresh,radius,microns)
		self.pca = self.chnls[primary_key].pcamed

		self.chnls[primary_key].align_data(self.chnls[primary_key].df_thresh,
			self.pca,comp_order,fit_dim,deg=deg)
		self.mm = self.chnls[primary_key].mm
		self.vertex = self.chnls[primary_key].vertex

		self.chnls[primary_key].transform_coordinates()

		print('Primary channel',primary_key,'processing complete')

		for ch in self.chnls.keys():
			if ch != primary_key:
				self.chnls[ch].preprocess_data(gthresh,scale,microns)
				
				self.chnls[ch].align_data(self.chnls[ch].df_thresh,
					self.pca,comp_order,fit_dim,deg=deg,
					mm = self.mm, vertex = self.vertex)

				self.chnls[ch].transform_coordinates()
				print(ch,'processed')

	def save_projections(self,subset):
		'''Save projections of both channels into files'''

		for ch in self.chnls.keys():
			fig = self.chnls[ch].plot_projections(self.chnls[ch].df_align,subset)
			fig.savefig(os.path.join(self.outdir,
				self.name+'_'+self.number+'_'+ch+'_MIP.png'))

		print('Projections generated')

	def save_psi(self):
		'''Save all channels into psi files'''

		columns = ['x','y','z','ac','r','theta']

		for ch in self.chnls.keys():
			write_data(os.path.join(self.outdir,
				self.name+'_'+self.number+'_'+ch+'.psi'),
				self.chnls[ch].df_align[columns])

		print('PSIs generated')

	def add_psi_data(self,filepath,key):
		''''Read in previously processed psi data as a dataframe'''

		self.chnls[key] = read_psi(filepath)

class math_model:
	'''Class to contain attributes and data associated with math model'''

	def __init__(self,model):

		self.cf = model
		self.p = np.poly1d(model)

def process_sample(num,root,outdir,name,chs,prefixes,threshold,scale,deg,primary_key,comp_order,fit_dim,flip_dim):
	'''Process single sample through embryo class and 
	save df to csv'''

	tic = time.time()

	print(num,'Starting',name)

	e = embryo(name,num,outdir)
	e.add_channel(os.path.join(root,chs[0],prefixes[0]+'_'+num+'_Probabilities.h5'),
		'at')
	e.add_channel(os.path.join(root,chs[1],prefixes[1]+'_'+num+'_Probabilities.h5'),'zrf1')
	e.process_channels(threshold,scale,deg,primary_key,[0,2,1],['x','z'],'z')
	
	print(num,'Data processing complete',name)

	e.save_projections(0.1)
	e.save_psi()

	toc = time.time()
	print(name,'complete',toc-tic)

##### PSI file processing ############

def write_header(f):
	'''Writes header for PSI file with columns Id,x,y,z,ac,r,theta'''

	contents = [
		'PSI Format 1.0',
		'',
		'column[0] = "Id"',
		'column[1] = "x"',
		'column[2] = "y"',
		'column[3] = "z"',
		'column[4] = "ac"',
		'symbol[4] = "A"',
		'type[4] = float',
		'column[5] = "r"',
		'symbol[5] = "R"',
		'type[5] = float',
		'column[6] = "theta"',
		'symbol[6] = "T"',
		'type[6] = float'
	]

	for line in contents:
		f.write('# '+ line + '\n')

def write_data(filepath,df):
	'''Writes data in PSI format to file after writing header'''

	#Open new file at given filepath
	f = open(filepath,'w')

	#Write header contents to file
	write_header(f)

	n = df.count()['x']
    
	#Write line with sample number
	f.write(str(n)+' 0 0\n')

	#Write translation matrix
	f.write('1 0 0\n'+
			'0 1 0\n'+
			'0 0 1\n')

	#Write dataframe to file using pandas to_csv function to format
	f.write(df.to_csv(sep=' ', index=True, header=False))

	f.close()

	print('Write to',filepath,'complete')

def read_psi(filepath):
	'''Reads psi file and saves data into dataframe'''

	df = pd.read_csv(filepath,
		sep=' ',
		header=19, #This value is now wrong
		names=['x','y','z','ac','r','theta']) #may also be wrong

	return(df)

def calculate_models(Ldf):
	'''Calculate models for each sample based on a list of at dataframes that are already aligned'''

	modeldf = pd.DataFrame({'a':[],'b':[],'c':[]})

	for df in Ldf:
		s = cranium.brain()
		s.add_aligned_df(df)
		s.fit_model(s.df_align,2)

		modeldf = modeldf.append(pd.DataFrame({'a':[s.mm.cf[0]],'b':[s.mm.cf[1]],'c':[s.mm.cf[2]]}))

	return(modeldf)

def fit_distribution(y,var,distname):
	'''Fit beta distribution to subset of data and return figure'''

	dist = getattr(scipy.stats, distname)
	sbst = y[var].sample(frac=0.01)

	x = np.linspace(y[var].min(),y[var].max(),100)

	param = dist.fit(sbst)
	if distname == 'beta':
		pdf_fit = dist.pdf(x,param[0],param[1],loc=param[2],scale=param[3])
	elif distname == 'uniform':
		pdf_fit = dist.pdf(x,param[0],param[1])

	fig,ax = plt.subplots()
	ax.hist(sbst,bins=15,normed=True)
	ax.plot(x,pdf_fit,c='y')
	#ax.set_ylim([0,1])

	ax.set_title(distname+str(param))

	return(fig,param)

def calculate_rms(df):
	'''Calculated the residual mean error for a sample'''

	return(np.sqrt(np.sum(df.r**2)/df.size))

def test_distributions(y,snum,plot=False,outdir=None):
	'''Data should be a single series already processed, e.g. sampled or absolute value'''

	#fit gamma distribution
	gout = fit_gamma(y,plot=plot)

	#fit beta distribuiton
	bout = fit_beta(y,plot=plot)

	#fit half normal distribution
	hnout = fit_norm(y,plot=plot)

	if plot==True:
		if outdir == None:
			outdir = os.getcwd()
		gout[-1].savefig(os.path.join(outdir,'gamma'+str(snum)+'.jpg'))
		bout[-1].savefig(os.path.join(outdir,'beta'+str(snum)+'.jpg'))
		hnout[-1].savefig(os.path.join(outdir,'halfnorm'+str(snum)+'.jpg'))

	df = pd.DataFrame({'snum':[snum],
		'g_a':[gout[0][0]], 'g_loc':[gout[0][1]], 'g_scale':[gout[0][2]],
		'g_D':[gout[1]],'g_p':[gout[2]],
		'b_a':[bout[0][0]], 'b_b':[bout[0][1]], 'b_loc':[bout[0][2]], 'b_scale':[bout[0][3]],
		'b_D':[bout[1]], 'b_p':[bout[2]],
		'hn_loc':[hnout[0][0]], 'hn_scale':[hnout[0][1]],
		'hn_D':[hnout[1]], 'hn_p':[hnout[2]]
		})

	return(df)

def test_beta(y,snum,plot=False,outdir=None):
	'''Data should be a single series already processed, e.g. sampled or absolute value'''

	#fit beta distribuiton
	bout = fit_beta(y,plot=plot)

	if plot==True:
		if outdir == None:
			outdir = os.getcwd()
		bout[-1].savefig(os.path.join(outdir,'beta'+str(snum)+'.jpg'))

	df = pd.DataFrame({'snum':[snum],
		'b_a':[bout[0][0]], 'b_b':[bout[0][1]], 'b_loc':[bout[0][2]], 'b_scale':[bout[0][3]],
		'b_D':[bout[1]], 'b_p':[bout[2]]
		})

	return(df)

def test_gamma(y,snum,plot=False,outdir=None):
	'''Data should be a single series already processed, e.g. sampled or absolute value'''

	#fit beta distribuiton
	gout = fit_gamma(y,plot=plot)

	if plot==True:
		if outdir == None:
			outdir = os.getcwd()
		gout[-1].savefig(os.path.join(outdir,'beta'+str(snum)+'.jpg'))

	df = pd.DataFrame({'snum':[snum],
		'g_a':[gout[0][0]], 'g_loc':[gout[0][1]], 'g_scale':[gout[0][2]],
		'g_D':[gout[1]],'g_p':[gout[2]]
		})

	return(df)
    
def fit_gamma(y,plot=False):
	'''Fit a gamma distribution to data
	Returns: parameters, ks test, and plot if specified'''

	gamma = scipy.stats.gamma

	#Calculate parameters of gamma distribution
	param = gamma.fit(y)

	#Calculate range over which to test data
	x = np.linspace(gamma.ppf(0.01,param[0],loc=param[1],scale=param[2]),
		gamma.ppf(0.99,param[0],loc=param[1],scale=param[2]),y.size)

	#Calculate pdf and cdf
	pdf = gamma.pdf(x,param[0],loc=param[1],scale=param[2])
	cdf = gamma.cdf(x,param[0],loc=param[1],scale=param[2])

	#Calculate estimated cdf
	ecdf_x = np.sort(y)
	ecdf_y = np.array(range(len(ecdf_x)))/float(len(ecdf_x))

	#Calculate ks test for goodness of fit
	D,pvalue = scipy.stats.kstest(y,'gamma',args=(param[0],param[1],param[2]))

	if plot == True:
		fig = plt.figure()
		ax = fig.add_subplot(121)
		ay = fig.add_subplot(122)

		ax.plot(ecdf_x,ecdf_y,label='ecdf')
		ax.plot(x,cdf,label='cdf')
		ax.set_title('Gamma distribution')
		ax.legend()

		ay.hist(y,normed=True)
		ay.plot(x,pdf)

		return(param,D,pvalue,fig)
	else:
		return(param,D,pvalue)

def fit_beta(y,plot=False):
	'''Fit beta distribution and test fit'''

	beta = scipy.stats.beta

	#Calculate beta parameters
	param = beta.fit(y)

	#Calculate x range
	x = np.linspace(beta.ppf(0.01,param[0],param[1],loc=param[2],scale=param[3]),
		beta.ppf(0.99,param[0],param[1],loc=param[2],scale=param[3]),y.size)

	#Calculate pdf and cdf
	pdf = beta.pdf(x,param[0],param[1],loc=param[2],scale=param[3])
	cdf = beta.cdf(x,param[0],param[1],loc=param[2],scale=param[3])

	#Calculate ecdf
	ecdf_x = np.sort(y)
	ecdf_y = np.array(range(len(ecdf_x)))/float(len(ecdf_x))

	#Calculate ks test for goodness of fit
	D,pvalue = scipy.stats.kstest(y,'beta',args=(param[0],param[1],param[2],param[3]))

	if plot == True:
		fig = plt.figure()
		ax = fig.add_subplot(121)
		ay = fig.add_subplot(122)

		ax.plot(ecdf_x,ecdf_y,label='ecdf')
		ax.plot(x,cdf,label='cdf')
		ax.set_title('Beta Distribution')
		ax.legend()

		ay.hist(y,normed=True)
		ay.plot(x,pdf)

		return(param,D,pvalue,fig)
	else:
		return(param,D,pvalue)

def fit_norm(y,plot=False):
	'''Fit halfnorm distribution'''

	halfnorm = scipy.stats.norm

	#Fit parameters
	param = halfnorm.fit(y)

	#Calculate x 
	x = np.linspace(halfnorm.ppf(0.01,loc=param[0],scale=param[1]),
		halfnorm.ppf(0.99,loc=param[0],scale=param[1]),y.size)

	#Calculate cdf and pdf
	pdf = halfnorm.pdf(x,loc=param[0],scale=param[1])
	cdf = halfnorm.cdf(x,loc=param[0],scale=param[1])

	#Calculate ecdf
	ecdf_x = np.sort(y)
	ecdf_y = np.array(range(len(ecdf_x)))/float(len(ecdf_x))

	#Calculate ks test for goodness of fit
	D,pvalue = scipy.stats.kstest(y,'halfnorm',args=(param[0],param[1]))

	if plot == True:
		fig = plt.figure()
		ax = fig.add_subplot(121)
		ay = fig.add_subplot(122)

		ax.plot(ecdf_x,ecdf_y,label='ecdf')
		ax.plot(x,cdf,label='cdf')
		ax.set_title('Halfnorm Fit')
		ax.legend()

		ay.hist(y,normed=True)
		ay.plot(x,pdf)

		return(param,D,pvalue,fig)
	else:
		return(param,D,pvalue)

def calculate_sample_error(y,param,distname,n):
	'''Calculate RMSE based on KDE of sample and PDF of distribution'''

	#Define x range
	x = np.linspace(np.min(y),np.max(y),n)

	#Calculate kde
	gkde = scipy.stats.gaussian_kde(y)
	kdepdf = gkde.evaluate(x)

	#Calculate pdf
	if distname == 'beta':
		pdf = scipy.stats.beta.pdf(x,param[0],param[1],loc=param[2],scale=param[3])
	elif distname == 'gamma':
		pdf = scipy.stats.gamma.pdf(x,param[0],loc=param[1],scale=param[2])

	error = np.sqrt(mean_squared_error(kdepdf,pdf))

	return(error)

def beta_pdf(x,param):
	return(scipy.stats.beta.pdf(x,param[0],param[1],loc=param[2],scale=param[3]))

def norm_pdf(x,param):
	return(scipy.stats.norm.pdf(x,loc=param[0],scale=param[1]))

def gamma_pdf(x,param):
	return(scipy.stats.gamma.pdf(x,param[0],loc=param[1],scale=param[2]))

def beta_rvs(param):
	return(scipy.stats.beta.rvs(param[0],param[1],loc=param[2],scale=param[3],size=50000))

def gamma_rvs(param):
	return(scipy.stats.gamma.rvs(param[0],loc=param[1],scale=param[2],size=50000))

def calculate_xcurve(ac,a):

	b=0

	xone = (-b + np.sqrt(b**2 + 4*a*ac))/2*a
	xtwo = (-b - np.sqrt(b**2 + 4*a*ac))/2*a

	if ac < 0:
		if xone < 0:
			return(xone)
		else:
			return(xtwo)
	else:
		if xone > 0:
			return(xone)
		else:
			return(xtwo)

def calculate_zcurve(x,a):

	p = np.poly1d([a,0,0])

	return(p(x))

def calculate_y(theta,r):

	return(r*np.sin(theta))

def calculate_xyz(row,a):

	ac,theta,r = row.ac,row.theta,row.r

	xc = calculate_xcurve(ac,a)
	zc = calculate_zcurve(xc,a)

	gamma = np.arctan(-1/(2*a*xc))
	k = r*np.cos(theta)
	alpha = k*np.sin(gamma)
	beta = k*np.cos(gamma)

	x = xc + alpha
	z = zc + beta
	y = calculate_y(theta,r)

	return(pd.Series({'x':x, 'y':y, 'z':z}))
	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

def generate_poc(ac,theta,r,a,outpath):
	'''dist name as last item in param list
	second to last is tuple with range of data'''

	if ac[-1] == 'beta':
		Lac = beta_rvs(ac[:-1])

	if theta[-1] == 'beta':
		Ltheta = beta_rvs(theta[:-1])

	if r[-1] == 'beta':
		Lr = beta_rvs(r[:-1])
	elif r[-1] == 'gamma':
		Lr = gamma_rvs(r[:-1])

	df = pd.DataFrame({
		'ac':Lac,
		'theta':Ltheta,
		'r':Lr
		})

	df = df.join(df.apply((lambda row: calculate_xyz(row,a)), axis=1))

	write_data(outpath,df)

def read_psi_to_dict(directory,dtype):
	'''Read psi in directory into dict of df'''

	dfs = {}
	for f in os.listdir(directory):
		if dtype in f:
			df = read_psi(os.path.join(directory,f))
			num = f.split('_')[-2][:2]
			dfs[num] = df

	return(dfs)

def concatenate_dfs(dfdict):
	'''Concatenated dfs from dict into one df'''

	L = []
	for key in dfdict.keys():
		L.append(dfdict[key])

	alldf = pd.concat(L)

	return(alldf)

def generate_kde(data,var,x,absv=False):
	'''Generate list of kde from either dict or list'''

	L = []

	#Dictionary workflow
	if type(data) ==  type({}):
		for key in data.keys():
			if absv == False:
				y = data[key][var]
			elif absv == True:
				y = np.abs(data[key][var])
			kde = scipy.stats.gaussian_kde(y).evaluate(x)
			L.append(kde)
	#List workflow
	elif type(data) == type([]):
		for df in data:
			if absv == False:
				y = df[var]
			elif absv == True:
				y = np.abs(df[var])
			kde = scipy.stats.gaussian_kde(y).evaluate(x)
			L.append(kde)

	return(L)

def fit_bimodal_theta(D,split,frac,x):
	'''Fit bimodal theta to concatenated df from dict D'''

	alldf = concatenate_dfs(D)
	y = alldf.theta.sample(frac=frac)

	lpoints = len(y[y<split])
	rpoints = len(y[y>split])

	pleft = fit_norm(y[y<split])
	pright = fit_norm(y[y>split])

	pdfl = norm_pdf(x,pleft[0])*(lpoints/(lpoints+rpoints))
	pdfr = norm_pdf(x,pright[0])*(rpoints/(lpoints+rpoints))

	return(pdfl+pdfr)

def calculate_area_error(pdf,Lkde,x):
	'''Calculate area between pdf and each Lkde'''

	L = []

	for kde in Lkde:
		L.append(simps(np.abs(pdf-kde),x))

	return(L)

def rescale_variable(Ddfs,var,newvar):
	'''Rescale variable from -1 to 1 and save in newvar in dict'''

	Dout = {}

	for key in Ddfs.keys():
		df = Ddfs[key]
		y = df[var]

		#Find min and max values
		mi,ma = y.min(),y.max()

		#Divide dataframe
		dfmin = df[y < 0]
		dfmax = df[y >= 0]

		#Add rescaled variable to dataframe
		dfmin = dfmin.join(pd.DataFrame({newvar:y/np.abs(mi)}))
		dfmax = dfmax.join(pd.DataFrame({newvar:y/np.abs(ma)}))

		dfout = dfmin.append(dfmax)
		Dout[key] = dfout

	return(Dout)