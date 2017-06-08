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

class brain:

	def __init__(self):
		'''Initialize brain object'''

	def read_data(self,filepath):
		'''Reads 3D data from file and selects appropriate channel'''

		#Read h5 file and extract probability data
		f = h5py.File(filepath,'r')
		c1 = np.array(f.get('exported_data')[:,:,:,0])
		c2 = np.array(f.get('exported_data')[:,:,:,1])

		#Figure out which channels has more zeros and therefore is background
		if np.count_nonzero(c1<0.1) > np.count_nonzero(c1>0.9):
			self.raw_data = c1
		else:
			self.raw_data = c2

	def create_dataframe(self):
		'''Creates a pandas dataframe containing the x,y,z and signal/probability value 
		for each point in the :py:attr:`brain.raw_data` array'''

		dim = self.raw_data.shape
		xyz = np.zeros((dim[0],dim[1],dim[2],4))

		#Generate array with xyz values for each point
		for x in range(dim[2]):

			xyz[:,:,x,2] = x
			xyz[:,:,x,3] = self.raw_data[:,:,x]

			zindex = np.arange(0,dim[0],1)
			yindex = np.arange(0,dim[1],1)
			gy,gz = np.meshgrid(yindex,zindex)

			xyz[:,:,x,0] = gz
			xyz[:,:,x,1] = gy

		flat = np.reshape(xyz,(-1,4))

		#Create dataframe of points
		self.df = pd.DataFrame({'x':flat[:,2],'y':flat[:,1],'z':flat[:,0],'value':flat[:,3]})

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

	def preprocess_data(self,threshold,scale):
		'''Thresholds and scales data prior to PCA'''

		#Create new dataframe with values above threshold
		self.threshold = threshold
		self.df_thresh = self.df[self.df.value > self.threshold]

		#Scale xyz by value in scale array to force PCA axis selection
		self.scale = scale
		self.df_scl = pd.DataFrame({
			'x':self.df_thresh.x * self.scale[0], 
			'y':self.df_thresh.y * self.scale[1],
			'z':self.df_thresh.z * self.scale[2]})

	def calculate_pca(self):
		'''Create pca object and calculate transformation matrix'''

		self.pca = PCA(n_components=3)
		self.pca.fit(self.df_scl[['x','y','z']])

	def add_pca(self,pca):
		'''Add pca object from another channel'''

		self.pca = pca

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
			self.fit_model(self.df_align,deg,fit_dim)
		else:
			self.mm = mm

	def fit_model(self,df,deg,fit_dim):
		'''Fit model to dataframe'''

		self.mm = math_model(np.polyfit(df[fit_dim[0]], df[fit_dim[1]], deg=deg))

	###### Functions associated with alpha, r, theta coordinate system ######

	def find_distance(self,t,point):
		'''Find euclidean distance between math model(t) and data point in the xy plane'''

		x = float(t)
		y = self.mm.p(x)

		#Calculate distance between two points passed as array
		dist = np.linalg.norm(point - np.array([x,y]))

		return(dist)

	def find_min_distance(self,row):
		'''Find the point on the curve that produces the minimum distance between the point and the data point'''

		dpoint = np.array([row.x,row.y])

		#Use scipy.optimize.minimize to find minimum solution of brain.find_distance
		result = minimize(self.find_distance, dpoint[0], args=(dpoint))

		x = result['x'][0]
		y = self.mm.p(x)
		z = 0
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

	def find_theta(self,row,xc,zc):
		'''Find theta value for a row describing angle between point and plane'''

		theta = np.arctan2(row.z-zc,row.x-xc)
		return(theta)

	def calc_coord(self,row):
		'''Calculate alpah, r, theta for a particular row'''

		xc,yc,zc,r = self.find_min_distance(row)
		ac = self.find_arclength(xc)
		theta = self.find_theta(row,xc,zc)

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

	def process_channels(self,threshold,scale,deg,primary_key,comp_order,fit_dim,flip_dim):
		'''Process channels through alignment'''

		#Process primary channel
		self.chnls[primary_key].create_dataframe()
		self.chnls[primary_key].preprocess_data(threshold,scale)
		self.chnls[primary_key].calculate_pca()
		self.chnls[primary_key].pca_transform(comp_order,fit_dim,flip_dim,deg=deg)
		self.chnls[primary_key].transform_coordinates()

		print('Primary channel',primary_key,'processing complete')

		for ch in self.chnls.keys():
			if ch != primary_key:
				self.chnls[ch].create_dataframe()
				self.chnls[ch].preprocess_data(threshold,scale)
				self.chnls[ch].add_pca(self.chnls[primary_key].pca)
				self.chnls[ch].pca_transform(comp_order,fit_dim,flip_dim,mm=self.chnls[primary_key].mm,
					flip=self.chnls[primary_key].flip,
					vertex=self.chnls[primary_key].vertex)
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