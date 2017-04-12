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

	def show_plane(dimension,plane):
		'''Shows specified plane'''

		if dimension=='x':
			data = sample.raw_data[:,:,plane]
		elif dimension=='y':
			data = sample.raw_data[:,plane,:]
		elif dimension=='z':
			data = sample.raw_data[plane,:,:]
		else:
			print('Invalid dimension specified, try "x","y","z"')

		fig,ax = plt.subplots()
		cax = ax.imshow(sample.raw_data[:,:,300],cmap='plasma')
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

		#Identify flat plane
		self.p_flat = smf.ols(formula='y ~ x + z',data=self.df_thresh).fit()

		#Identify parabolic plane
		self.p_para = smf.ols(formula='z ~ y + x + I(x**2)',data=self.df_thresh).fit()

		#Find intersection
		self.model = {
		'a' : self.p_para.params[3],
		'b' : self.p_para.params[2],
		'c' : self.p_para.params[1],
		'd' : self.p_para.params[0],
		'e' : self.p_flat.params[1],
		'f' : self.p_flat.params[2],
		'g' : self.p_flat.params[0]
		}

		self.model['a_prime'] = (self.model['a']*self.model['f']) / (1 - self.model['c']*self.model['f'])
		self.model['b_prime'] = (self.model['e'] + self.model['b']*self.model['f']) / (1 - self.model['c']*self.model['f'])
		self.model['c_prime'] = (self.model['g'] + self.model['d']*self.model['f']) / (1 - self.model['c']*self.model['f'])