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

	###### Functions associated with alpha, r, theta coordinate system ######

	def find_distance(self,t,point):
		'''Find euclidean distance between point on line defined by t and data point'''

		x = float(t)
		y = self.mm.calc_y(x)
		z = self.mm.calc_z(x)

		dist = np.linalg.norm(point - np.array([x,y,z]))

		return(dist)

	def find_min_distance(self,row):
		'''Find the point on the curve that produces the minimum distance between the point and the data point'''

		dpoint = np.array([row.x,row.y,row.z])

		#Use scipy.optimize.minimize to find minimum solution of brain.find_distance
		result = minimize(self.find_distance, dpoint[0],args=(dpoint))

		x = result['x'][0]
		y = self.mm.calc_y(x)
		z = self.mm.calc_z(x)
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
	'''Class to contain attribues and data associated with math model'''

	def __init__(self,coef,p,x,y,z):

		self.coef = coef
		self.p = p
		self.x = x
		self.y = y
		self.z = z

		self.find_vertex()
		self.find_focus()

	def calc_y(self,t):
		'''Calculate y value according to a given t'''

		y = self.p['ay']*(t**2) + self.p['by']*t + self.p['cy']
		return(y)

	def calc_z(self,t):
		'''Calculate z value according to a given t'''

		z = self.p['az']*(t**2) + self.p['bz']*t + self.p['cz']
		return(z)

	def find_vertex(self):
		'''Calculate position of vertex'''

		self.vx = -self.coef['b']/(2*self.coef['a'])
		self.vy = self.calc_y(self.vx)
		self.vz = self.calc_z(self.vx)

	def find_focus(self):
		'''Calculate position of focus'''

		dy = np.sqrt(1/(16*self.coef['a']*((1/self.coef['f'])**2 + 1)))
		dz = (1/self.coef['f'])*dy

		self.fx = self.vx
		self.fy = self.vy + dy
		self.fz = self.vz + dz

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