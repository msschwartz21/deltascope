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
from sklearn.preprocessing import scale
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

	def plot_projections(self,df):
		'''Plots x, y, and z projections from the dataframe'''
    
    	#Create figure and subplots
		fig = plt.figure(figsize=(12,6))
		ax = fig.add_subplot(131)
		ay = fig.add_subplot(132)
		az = fig.add_subplot(133)

		#Create scatter plot for each projection
		ax.scatter(df['x'],df['z'])
		ay.scatter(df['x'],df['y'])
		az.scatter(df['z'],df['y'])

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
		plt.show()

	def align_sample(self,threshold,scale,deg):
		'''Realigns sample axes using PCA and translates so that the vertex is at the origin'''

		#Create new dataframe with values above threshold
		self.threshold = threshold
		self.df_thresh = self.df[self.df.value > self.threshold]

		#Scale xyz by value in scale array to force PCA axis selection
		self.df_scl = pd.DataFrame({
			'x':self.df_thresh.x * scale[0], 
			'y':self.df_thresh.y * scale[1],
			'z':self.df_thresh.z * scale[2]})

		#Fit pca to data and transform data points
		pca = PCA(n_components=3)
		pca_fit = pca.fit_transform(self.df_scl[['x','y','z']])

		#Create pca dataframe and remove scaling
		df_unscl = pd.DataFrame({
			'x':pca_fit[:,0]/scale[0],
			'y':pca_fit[:,1]/scale[1],
			'z':pca_fit[:,2]/scale[2]
			})

		#Find first model
		model = np.polyfit(df_unscl['x'],df_unscl['y'],deg=deg)
		p = np.poly1d(model)

		#Find vertex
		vx = -model[1]/(2*model[0])
		vy = p(vx)
		vz = df_unscl.z.mean()

		#Translate data so that the vertex is at the origin
		self.df_align = pd.DataFrame({
			'x': df_unscl.x - vx,
			'y': df_unscl.y - vy,
			'z': df_unscl.z - vz
			})

		#Calculate final model based on translated and aligned data
		self.model = math_model(np.polyfit(self.df_align.x, self.df_align.y,deg=deg))

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
		y = self.mm.pp(x)
		z = 0
		r = result['fun']

		return(x,y,z,r)
		#return(pd.Series({'xc':x, 'yc':y, 'zc':z, 'r':r}))

	def find_alpha(self,xc,yc,zc):
		'''Calculate alpha for a row containing point data'''

		#Calculate distances between vertex, focus, and point on the curve
		vf = np.linalg.norm(np.array([self.mm.vx,self.mm.vy,self.mm.vz]) - 
			np.array([self.mm.fx,self.mm.fy,self.mm.fz]))
		fp = np.linalg.norm(np.array([self.mm.fx,self.mm.fy,self.mm.fz]) - 
			np.array([xc,yc,zc]))
		pv = np.linalg.norm(np.array([xc,yc,zc]) - 
			np.array([self.mm.vx,self.mm.vy,self.mm.vz]))

		alpha = np.arccos((vf**2 + fp**2 - pv**2)/(2*vf*fp))

		#Set alpha sign based on position along x axis in relationship to vertex
		if xc >= self.mm.vx:
			return(alpha)
		elif xc < self.mm.vx:
			return(-alpha)

	def integrand(self,x):
		'''Function to integrate to calculate arclength'''

		y_prime = self.mm.p['ay']*2*x + self.mm.p['by']
		z_prime = self.mm.p['az']*2*x + self.mm.p['bz']

		arclength = np.sqrt(y_prime**2 + z_prime**2)
		return(arclength)

	def find_length(self,xc):
		'''Calculate arclength for a row'''

		ac,err = scipy.integrate.quad(self.integrand,xc,self.mm.vx)
		return(ac)

	def dist_to_plane(self,xz,row):
		'''Find shortest distance between point and the plane'''

		#Find y value based on a given x and z
		y = self.mm.coef['f']*xz[1] + self.mm.coef['e']*xz[0] + self.mm.coef['g']

		dist = np.linalg.norm(np.array([xz[0],y,xz[1]]) - 
			np.array([row.x,row.y,row.z]))
		return(dist)

	def find_theta(self,row,r,zc):
		'''Find theta value for a row describing angle between point and plane'''

		#Find the point which minimizes the distance between the point and the plane
		result = minimize(self.dist_to_plane,[row.x,row.z],args=(row))
		planey = self.mm.coef['f']*result['x'][1] + self.mm.coef['e']*result['x'][0] + self.mm.coef['g']

		#Set sign of r based on whether its left or right of y axis
		if row.z >= zc:
			r = r
		elif row.z < zc:
			r = -r

		theta = np.arccos(result['fun']/r)

		#Change sign of theta based on whether its above or below the plane
		ud = self.mm.coef['e']*row.x + self.mm.coef['f']*row.z - row.y + self.mm.coef['g']
		if ud >= 0: 
			return(theta)
		elif ud < 0:
			return(-theta)

	def calc_coord(self,row):
		'''Calculate alpah, r, theta for a particular row'''

		xc,yc,zc,r = self.find_min_distance(row)
		ac = self.find_length(xc)
		theta = self.find_theta(row,r,zc)

		return(pd.Series({'xc':xc, 'yc':yc, 'zc':zc,
					'r':r, 'ac':ac, 'theta':theta}))

	def transform_coordinates(self):
		'''Transform coordinate system so that each point is defined relative to 
		math model by (alpha,theta,r) (only applied to df_thresh'''

		#Calculate alpha, theta, r for each row in dataset
		self.df_thresh = self.df_thresh.join(self.df_thresh.apply((lambda row: self.calc_coord(row)), axis=1))

	def subset_data(self,sample_frac=0.5):
		'''Subset data based on proportion set in sample_frac'''

		self.subset = self.df_thresh.sample(frac=sample_frac)

	def add_thresh_df(self,df):
		'''Add dataframe of thresholded and transformed data to self.df_thresh'''

		self.df_thresh = df

class math_model:
	'''Class to contain attributes and data associated with math model'''

	def __init__(self,model):

		self.cf = model
		self.p = np.poly1d(model)

def process_sample(filepath):
	'''Process single sample through brain class and 
	save df to csv'''

	tic = time.time()

	path = '\\'.join(filepath.split('\\')[:-1])
	name = filepath.split('\\')[-1].split('.')[0]

	print('Starting',name)

	s = brain()
	s.read_data(filepath)
	s.create_dataframe()
	s.fit_model(0.9)
	s.transform_coordinates()
	s.df_thresh.to_csv(os.path.join(path,name+'.csv'))

	toc = time.time()
	print(name,'complete',toc-tic)