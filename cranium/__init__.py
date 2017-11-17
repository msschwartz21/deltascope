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
			c1 = np.array(d[:,:,:,9])
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

		#NB: scale variable actually contains microns dimensions

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

	def calculate_pca_median_2d(self,data,threshold,radius,microns):

		self.median = self.process_alignment_data(data,threshold,radius,microns)

		self.pcamed = PCA()
		self.pcamed.fit(self.median[['y','z']])

	def pca_transform_2d(self,df,pca,comp_order,fit_dim,deg=2,mm=None,vertex=None,flip=None):

		fit = pca.transform(df[['y','z']])
		df_fit = pd.DataFrame({
			'x':df.x,
			'y':fit[:,comp_order[1]-1],
			'z':fit[:,comp_order[2]-1]
			})

		self.align_data(df_fit,fit_dim,deg=2,mm=None,vertex=None,flip=None)

	def pca_transform_3d(self,df,pca,comp_order,fit_dim,deg=2,mm=None,vertex=None,flip=None):

		fit = pca.transform(df[['x','y','z']])
		df_fit = pd.DataFrame({
			'x':fit[:,comp_order[0]],
			'y':fit[:,comp_order[1]],
			'z':fit[:,comp_order[2]]
			})

		self.align_data(df_fit,fit_dim,deg=2,mm=None,vertex=None,flip=None)
	
	def align_data(self,df_fit,fit_dim,deg=2,mm=None,vertex=None,flip=None):
		'''Apply PCA transformation matrix and align data so that 
		the vertex is at the origin'''
		
		#If vertex for translation is not included
		if vertex == None and deg==2:
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
		elif deg == 1:
			#Calculate model
			model = np.polyfit(df_fit[fit_dim[0]],df_fit[fit_dim[1]],deg=deg)
			p = np.poly1d(model)

			if fit_dim[0] == 'x':
				vx = df_fit.x.mean()
				if fit_dim[1] == 'y':
					vy = p(vx)
					vz = df_fit.z.mean()
				else:
					vz = p(vx)
					vy = df_fit.y.mean()
			elif fit_dim[0] == 'y':
				vy = df_fit.y.mean()
				if fit_dim[1] == 'x':
					vx = p(vy)
					vz = df_fit.z.mean()
				else:
					vz = p(vy)
					vx = df_fit.x.mean()
			elif fit_dim[0] == 'z':
				vz = df_fit.z.mean()
				if fit_dim[1] == 'x':
					vx = p(vz)
					vy = df_fit.y.mean()
				else:
					vy = p(vz)
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
		if flip == None or flip == False:
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
		#r = result['fun']

		return(x,y,z)
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
		theta = np.arctan2(row.y-yc,row.z-zc)
		return(theta)

	def find_r(self,row,zc,yc):

		r = np.sqrt((row.z-zc)**2 + (row.y-yc)**2)
		return(r)

	def calc_coord(self,row):
		'''Calculate alpah, r, theta for a particular row'''

		xc,yc,zc = self.find_min_distance(row)
		ac = self.find_arclength(xc)
		theta = self.find_theta(row,zc,yc)
		r = self.find_r(row,zc,yc)

		return(pd.Series({'x':row.x,'y':row.y,'z':row.z,'xc':xc, 'yc':yc, 'zc':zc,
					'r':r, 'ac':ac, 'theta':theta}))

	def transform_coordinates(self):
		'''Transform coordinate system so that each point is defined relative to 
		math model by (alpha,theta,r) (only applied to df_align'''

		#Calculate alpha, theta, r for each row in dataset
		self.df_align = self.df_align.merge(self.df_align.apply((lambda row: self.calc_coord(row)), axis=1))

	def subset_data(self,sample_frac=0.5):
		'''Subset data based on proportion set in sample_frac'''

		self.subset = self.df_align.sample(frac=sample_frac)

	def add_thresh_df(self,df):
		'''Add dataframe of thresholded and transformed data to self.df_thresh'''

		self.df_thresh = df

	def add_aligned_df(self,df):
		'''Add dataframe of aligned data'''

		self.df_align = df
		self.mm = self.fit_model(self.df_align,2,['x','z'])

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

class landmarks:

	def __init__(self,percbins=[10,50,90],rnull=15):

		self.lm_wt_rf = pd.DataFrame()
		self.lm_mt_rf = pd.DataFrame()

		self.rnull = rnull
		self.percbins = percbins

	def calc_bins(self,Ldf,ac_num,tstep):
		'''Calculates alpha and theta bins based on ac_num and tstep'''

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
			self.acbins = np.linspace(-acmax,acmax)

		#Calculate tbins divisions based on tstep
		self.tbins = np.arange(-np.pi,np.pi+tstep,tstep)

	def calc_perc(self,df,snum,dtype,out):
		'''Calculate landmarks for a dataframe based on the bins and percentiles that have been previously defined'''

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
	'''Take a dataframe in which columns contain the bin parameters and 
	convert to a cartesian coordinate system'''

    ndf = pd.DataFrame()
    for c in df.columns:
        if len(c.split('_')) == 6:
            amn,amx,tmn,tmx,p,dtype = c.split('_')
            x = np.mean([float(amn),float(amx)])
            t = np.mean([float(tmn),float(tmx)])
            
            if dtype == 'r':
                r = np.mean(df[c])
                r_std = np.std(df[c])
                y = np.sin(t)*r
                z = np.sin(t)*r
                
                pts = np.mean(df['_'.join([amn,amx,tmn,tmx,p,'pts'])])
                
                D = pd.Series({'x':x,'y':y,'z':z,'r':r,'r_std':r_std,'t':t,'pts':pts})
                
                ndf = ndf.append(pd.Series(D),ignore_index=True)
            
    return(ndf)

def convert_to_arr(xarr,tarr,wt,mt):
	wtarr = np.zeros((len(xarr),len(tarr),wt.count(axis=0)['Unnamed: 0']))
	mtarr = np.zeros((len(xarr),len(tarr),wt.count(axis=0)['Unnamed: 0']))

	for c in mt.columns:
		if len(c.split('_') == 6):
			amn,amx,tmn,tmx,p,dtype = c.split('_')
			x = np.mean([float(amn),float(amx)])
			t = np.mean([float(tmn),float(tmx)])

			if dtype=='r':
				wtarr[np.where(xarr==x)[0],np.where(tarr==t)[0]] = wt[c]
				mtarr[np.where(xarr==x)[0],np.where(tarr==t)[0]] = mt[c]

	return(wtarr,mtarr)

P = {
	'zln':2,'zpt':3,'zfb':1,
	'wtc':'b','mtc':'r',
	'alpha':0.3,
	'cmap':'Greys_r',
	'xarr':None,
	'tarr':None
}

def subplot_lmk(ax,p,avg,sem,parr,xarr,tarr,dtype,Pn=P):

	for k in Pn.keys():
		P[k] = Pn[k]

	if dtype=='wt':
		c = P['wtc']
	else:
		c = P['mtc']

	ti1 = np.where(tarr==p[0])[0][0]
	ti2 = np.where(tarr==p[1])[0][0]

	ax.fill_between(xarr,avg[:,ti1]+sem[:,ti1],avg[:,ti1]-sem[:,ti1],alpha=P['alpha'],color=c,zorder=P['zfb'])
	ax.fill_between(xarr,-avg[:,ti2]+sem[:,ti2],-avg[:,ti2]-sem[:,ti2],alpha=P['alpha'],color=c,zorder=P['zfb'])

	ax.plot(xarr,avg[:,ti1],c=c,zorder=P['zln'])
	ax.plot(xarr,-avg[:,ti2],c=c,zorder=P['zln'])

	if dtype=='mt':
		ax.scatter(xarr,avg[:,ti1],c=parr[:,ti1],cmap=P['cmap'],zorder=P['zpt'])
		ax.scatter(xarr,-avg[:,ti2],c=parr[:,ti2],cmap=P['cmap'],zorder=P['zpt'])

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
	try:
		f.write(df[['x','y','z','ac','theta','r']].to_csv(sep=' ', index=True, header=False))
	except:
		f.write(df[['x','y','z']].to_csv(sep=' ', index=True, header=False))

	f.close()

	print('Write to',filepath,'complete')

def read_psi(filepath):
	'''Reads psi file and saves data into dataframe'''

	df = pd.read_csv(filepath,
		sep=' ',
		header=19)

	if len(df.columns) <= 4:
		df.columns = ['i','x','y','z']
	elif len(df.columns) >= 6:
		df.columns=['i','x','y','z','ac','theta','r']

	return(df)

def read_psi_to_dict(directory,dtype):
	'''Read psi in directory into dict of df'''

	dfs = {}
	for f in os.listdir(directory):
		if dtype in f:
			df = read_psi(os.path.join(directory,f))
			num = f.split('_')[-1][:-4]
			dfs[num] = df

	return(dfs)

###### Stand alone functions

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

def calculate_models(Ldf):
	'''Calculate models for each sample based on a list of at dataframes that are already aligned'''

	modeldf = pd.DataFrame({'a':[],'b':[],'c':[]})

	for df in Ldf:
		s = cranium.brain()
		s.add_aligned_df(df)
		s.fit_model(s.df_align,2)

		modeldf = modeldf.append(pd.DataFrame({'a':[s.mm.cf[0]],'b':[s.mm.cf[1]],'c':[s.mm.cf[2]]}))

	return(modeldf)

def generate_kde(data,var,x,absv=False):
	'''Generate list of kde from either dict or list'''

	L = []

	#Dictionary workflow
	if type(data) ==  type({}):
		for key in data.keys():
			if absv == False:
				y = data[key][var].sample(frac=0.1)
			elif absv == True:
				y = np.abs(data[key][var].sample(frac=0.1))
			kde = scipy.stats.gaussian_kde(y).evaluate(x)
			L.append(kde)
	#List workflow
	elif type(data) == type([]):
		for df in data:
			if absv == False:
				y = df[var].sample(frac=0.1)
			elif absv == True:
				y = np.abs(df[var].sample(frac=0.1))
			kde = scipy.stats.gaussian_kde(y).evaluate(x)
			L.append(kde)

	return(L)

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