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
import plotly.plotly as py
import plotly.graph_objs as go
from scipy.optimize import minimize


class brain:

	def __init__(self,filepath):
		'''Read raw data'''

		self.raw_data = self.read_data(filepath)

	def read_data(self,filepath):
		'''Reads 3D data from file and selects appropriate channel'''

		#Read h5 file and extract probability data
		f = h5py.File(filepath,'r')
		c1 = np.array(f.get('exported_data')[:,:,:,0])
		c2 = np.array(f.get('exported_data')[:,:,:,1])

		#Figure out which channels has more zeros and therefore is background
		if np.count_nonzero(c1<0.1) > np.count_nonzero(c1>0.9):
			return(c1)
		else:
			return(c2)

	def show_plane(self,dimension,plane):
		'''Shows specified plane'''

		if dimension=='x':
			data = self.raw_data[:,:,plane]
		elif dimension=='y':
			data = self.raw_data[:,plane,:]
		elif dimension=='z':
			data = self.raw_data[plane,:,:]
		else:
			print('Invalid dimension specified, try "x","y","z"')

		fig,ax = plt.subplots()
		cax = ax.imshow(data,cmap='plasma')
		fig.colorbar(cax)

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

	def fit_model(self,threshold):
		'''Calculates the mathematical model of the data'''

		self.threshold = threshold
		self.df_thresh = self.df[self.df.value > self.threshold]

		#Create xyz arrays for range of data
		x = np.linspace(self.df_thresh.x.min(),self.df_thresh.x.max())
		y = np.linspace(self.df_thresh.y.min(),self.df_thresh.y.max())
		z = np.linspace(self.df_thresh.z.min(),self.df_thresh.z.max())

		#Identify flat plane
		flat_model = smf.ols(formula='y ~ x + z',data=self.df_thresh).fit()
		xx,zz = np.meshgrid(x,z)
		Y = flat_model.params[0] + flat_model.params[1]*xx + flat_model.params[2]*zz
		self.f_plane = plane(flat_model,xx,Y,zz)

		#Identify parabolic plane
		para_model = smf.ols(formula='z ~ y + x + I(x**2)',data=self.df_thresh).fit()
		xx,yy = np.meshgrid(x,y)
		Z = para_model.params[0] + para_model.params[1]*yy + para_model.params[2]*xx + para_model.params[3]*(xx**2)
		self.p_plane = plane(para_model,xx,yy,Z)

		#Find intersection
		# y = ex + fz + g
		# z = ax**2 + bx + cy + d
		model = {
		'a' : para_model.params[3],
		'b' : para_model.params[2],
		'c' : para_model.params[1],
		'd' : para_model.params[0],
		'e' : flat_model.params[1],
		'f' : flat_model.params[2],
		'g' : flat_model.params[0]
		}

		model['a_prime'] = (model['a']*model['f']) / (1 - model['c']*model['f'])
		model['b_prime'] = (model['e'] + model['b']*model['f']) / (1 - model['c']*model['f'])
		model['c_prime'] = (model['g'] + model['d']*model['f']) / (1 - model['c']*model['f'])

		p = {
		'ay' : model['a_prime'],
		'by' : model['b_prime'],
		'cy' : model['c_prime'],
		'az' : model['a'] + model['c']*model['a_prime'],
		'bz' : model['b'] + model['c']*model['b_prime'],
		'cz' : model['d'] + model['c']*model['c_prime']
		}

		#Parametric equation in terms of t
		t = np.arange(0,1000)
		x_line = t
		y_line = model['a_prime']*(t**2) + model['b_prime']*t + model['c_prime']
		z_line = ((model['a'] + model['c']*model['a_prime'])*(t**2) + 
					(model['b'] + model['c']*model['b_prime'])*t + 
					model['c']*model['c_prime'] + model['d'])

		self.mm = math_model(model,p,x_line,y_line,z_line)

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

		return(alpha)

	def dist_to_plane(self,xz,row):
		'''Find shortest distance between point and the plane'''

		#Find y value based on a given x and z
		y = self.mm.coef['f']*xz[1] + self.mm.coef['e']*xz[0] + self.mm.coef['g']

		dist = np.linalg.norm(np.array([xz[0],y,xz[1]]) - 
			np.array([row.x,row.y,row.z]))
		return(dist)

	def find_theta(self,row,r):
		'''Find theta value for a row describing angle between point and plane'''

		result = minimize(self.dist_to_plane,[row.x,row.z],args=(row))

		theta = np.arccos(result['fun']/r)
		return(theta)

	def calc_coord(self,row):
		'''Calculate alpah, r, theta for a particular row'''

		xc,yc,zc,r = self.find_min_distance(row)
		alpha = self.find_alpha(xc,yc,zc)
		theta = self.find_theta(row,r)

		return(pd.Series({'xc':xc, 'yc':yc, 'zc':zc,
					'r':r, 'alpha':alpha, 'theta':theta}))

	def transform_coordinates(self):
		'''Transform coordinate system so that each point is defined relative to 
		math model by (alpha,theta,r) (only applied to df_thresh'''

		#Calculate alpha, theta, r for each row in dataset
		self.df_thresh = self.df_thresh.join(self.df_thresh.apply((lambda row: self.calc_coord(row)), axis=1))

	def subset_data(self,sample_frac=0.5):
		'''Subset data based on proportion set in sample_frac'''

		self.subset = self.df_thresh.sample(frac=sample_frac)

	def plot_model(self,sample_frac=0.5):
		'''Plot two planes, line model, and percentage of points
		Returns plotly fig object to be displayed according to user preference'''

		subset = self.df_thresh.sample(frac=sample_frac)

		points = dict(
			x = subset.x,
			y = subset.y,
			z = subset.z,
			type = 'scatter3d',
			mode = 'markers',
			marker = dict(
				size=3,
				color='black',
				opacity=0.01
				)
			)

		line = dict(
			x = self.mm.x[:600],
			y = self.mm.y[:600],
			z = self.mm.z[:600],
			type = 'scatter3d',
			mode = 'lines',
			line = dict(
				width = 3,
				color = 'green'
				)
			)

		flat = dict(
			x = self.f_plane.xx,
			y = self.f_plane.yy,
			z = self.f_plane.zz,
			type = 'surface',
			opacity = 0.6,
			showscale = False
			)

		para = dict(
			x = self.p_plane.xx,
			y = self.p_plane.yy,
			z = self.p_plane.zz,
			type = 'surface',
			opacity = 0.6,
			showscale = False
			)

		data = [points,line,flat,para]

		layout = dict(
			margin = dict(l=0,r=0,b=0,t=0)
			)

		fig = go.Figure(data=data,layout=layout)

		return(fig)


class plane:
	'''Class to contain attributes and data associated with a plane'''

	def __init__(self,model,xx,yy,zz):
		'''Save global variables'''

		self.model = model
		self.xx = xx
		self.yy = yy
		self.zz = zz

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