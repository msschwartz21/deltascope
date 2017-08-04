#!/usr/bin/env python

import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

import os
import yaml
import copy
import cranium
import numpy as np
import time

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class channelMaster:
	def __init__(self,parent,text,n):
		self.name = tk.Entry(parent)
		self.name.grid(row=n,column=1)
		self.name.insert(0,text)
		self.dir = None
		self.dButton = ttk.Button(parent,text='Select Folder',command=self.open_dir)
		self.dButton.grid(row=n,column=2)
		self.files = None

	def open_dir(self):
		self.dir = filedialog.askdirectory()
		self.dButton['text'] = self.dir

	def return_data_dict(self):
		out = {}
		#out['files'] = self.files
		out['dir'] = self.dir
		out['name'] = self.name
		print(out)

		return(out)

class Sample:
	def __init__(self,key):
		self.key = key
		self.fs = ['','','','']
		self.outpaths = ['','','','']

		self.cbVar = tk.IntVar()
		self.degVar = tk.StringVar()
		self.directVar = tk.StringVar()

	def tab2row(self,tab,row):

		self.checkB = tk.Checkbutton(tab,text=self.key,variable=self.cbVar)
		self.checkB.select()
		self.checkB.grid(row=row,column=0)

		for i,f in enumerate(self.fs):
			ttk.Label(tab,text=f).grid(row=row,column=i+1)

	def tab3row(self,tab,r,app):

		#Sample number label
		tk.Label(tab,text=self.key).grid(row=r,column=0)

		self.plotB = tk.Button(tab,text='Plot',command=lambda:self.first_alignment(app.p,app.Lc))
		self.plotB.grid(row=r,column=1)

		rotDegrees = ('-180','-90','0','90','180')
		rotDirection = ('X','Y','Z')

		self.degVar.set(rotDegrees[2])
		self.degM = tk.OptionMenu(tab,self.degVar,*rotDegrees)
		self.degM.grid(row=r,column=2)

		self.directVar.set(rotDirection[0])
		self.directM = tk.OptionMenu(tab,self.directVar,*rotDirection)
		self.directM.grid(row=r,column=3)



	def first_alignment(self,p,Lc):
		tic = time.time()

		pc = self.check_settings(p)

		e = cranium.embryo(pc['expname'],self.key,pc['outdir'])
		for i,c in enumerate(Lc):
			if c.dir != None:
				e.add_channel(os.path.join(c.dir,self.fs[i]),c.name)
				e.chnls[c.name].preprocess_data(pc['genthresh'],[1,1,1],pc['microns'])

				#Processing for the structural channel
				if i==0:
					e.chnls[c.name].calculate_pca_median(e.chnls[c.name].raw_data,pc['medthresh'],
						pc['radius'],pc['microns'])
					pca = e.chnls[c.name].pcamed
					e.chnls[c.name].align_data(e.chnls[c.name].df_thresh,pca,
						pc['comporder'],pc['fitdim'],deg=pc['deg'])
					mm = e.chnls[c.name].mm
					vertex = e.chnls[c.name].vertex

				#Secondary channel processing
				e.chnls[c.name].align_data(e.chnls[c.name].df_thresh,pca,
					pc['comporder'],pc['fitdim'],
					deg=pc['deg'],mm=mm,vertex=vertex)

			print('sample alignment complete',time.time()-tic)

			#Save data
			for i,c in enumerate(Lc):
				if c.dir != None:
					cname = c.name.get()
					for p in ':,./\\[]{};() ':
						cname = cname.replace(p,'')

					fpath = os.path.join(pc['outdir'],'_'.join([cname,pc['expname'],self.key])+'.psi')
					cranium.write_data(fpath,e.chnls[c.name].df_align)
					self.outpaths[i] = fpath

			self.plot_projections([e.chnls[Lc[0].name].df_thresh, e.chnls[Lc[0].name].df_align],mm)

	def plot_projections(self,Ldf,mm):
		n = len(Ldf)
		fig, axarr = plt.subplots(n,3,num='Sample '+self.key) #(rows, columns)

		for i,df in enumerate(Ldf):
			df = df.sample(frac=0.01)
			axarr[i,0].scatter(df.x,df.z)
			axarr[i,1].scatter(df.x,df.y)
			axarr[i,2].scatter(df.z,df.y)

			##### To plot the model I need to account for which dimensions the model is supposed to be fit in 

			axarr[i,0].set_title('Y projection')
			axarr[i,1].set_title('Z projection')
			axarr[i,2].set_title('X projection')
			axarr[i,0].set_xlabel('X')
			axarr[i,0].set_ylabel('Z')
			axarr[i,1].set_xlabel('X')
			axarr[i,1].set_ylabel('Y')
			axarr[i,2].set_xlabel('Z')
			axarr[i,2].set_ylabel('Y')

		fig.show()

	def check_settings(self,p):

		pc = {}

		#Validate experiment name and outdir
		v = p['expname']['entry'].get()
		if v == '':
			messagebox.showerror('Error','Experiment name is not defined')
			return
		else:
			pc['expname'] = v

		if self.outdir == None:
			messagebox.showerror('Error','Output directory has not been selected')
			return
		else:
			pc['outdir'] = self.outdir

		#Check median threshold value
		s = p['medthresh']['entry'].get()
		try:
			v = float(s)
			if v <= 1 and v >= 0:
				pc['medthresh'] = v
			else:
				messagebox.showerror('Error','Median threshold input must be a number between 0 and 1')
		except ValueError:
			messagebox.showerror('Error','Median threshold input must be a number between 0 and 1')
			return

		#Check radius must be integer
		s = p['radius']['entry'].get()
		try:
			pc['radius'] = int(s)
		except ValueError:
			messagebox.showerror('Error','Radius input must be an integer')
			return

		#Check genthresh float between 0 and 1
		s = p['genthresh']['entry'].get()
		try:
			v = float(s)
			if v <= 1 and v >= 0:
				pc['genthresh'] = v
			else:
				messagebox.showerror('Error','General threshold input must be a number between 0 and 1')
				return
		except ValueError:
			messagebox.showerror('Error','General threshold input must be a number between 0 and 1')
			return

		#Check microns
		s = p['microns']['entry'].get()
		ls = s.replace(' ','').split(',')
		if len(ls) != 3:
			messagebox.showerror('Error','Micron input must be a list of three numeric values seperated by commas')
			return
		else:
			lv = []
			for m in ls:
				try:
					v = float(m)
					lv.append(v)
				except ValueError:
					messagebox.showerror('Error','Micron input must be a list of three numeric values seperated by commas')
					return
			pc['microns'] = lv

		#Check degree
		s = p['deg']['entry'].get()
		try:
			pc['deg'] = int(s)
		except ValueError:
			messagebox.showerror('Error','Degree input must be an integer')
			return

		#Check component order
		s = p['comporder']['entry'].get()
		ls = s.replace(' ','').split(',')
		if len(ls) != 3:
			messagebox.showerror('Error','Component order input must be a list of three integer values seperated by commas')
			return
		else:
			lv = []
			for m in ls:
				try:
					v = int(m)
					lv.append(v)
				except ValueError:
					messagebox.showerror('Error','Component input must be a list of three integer values seperated by commas')
					return
			pc['comporder'] = lv

		#Check fit dimensions
		s = p['fitdim']['entry'].get()
		ls = s.replace(' ','').split(',')
		if len(ls) != 2:
			messagebox.showerror('Error','Fit dimensions must be two of x,y or z seperated by a comma')
			return
		else:
			lv = []
			for m in ls:
				if m in ['x','y','z']:
					lv.append(m)
				else:
					messagebox.showerror('Error','Fit dimensions must be two of x,y or z seperated by a comma')
					return
			pc['fitdim'] = lv

		return(pc)
