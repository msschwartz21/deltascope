import cranium
import multiprocessing as mp
import time
import os
from functools import partial
from sys import argv

def check_config(D,path):

	#Validate that experiment name is defined
	if D['expname'] == '':
		print('Experiment name (expname) is not defined in',path)
		raise
	else:
		expname = D['expname']

	#Check median threshold value
	if D['medthresh'] <= 1 and D['medthresh'] >= 0:
		medthresh = D['medthresh']
	else:
		print('Median threshold (medthresh) must be a number between 0 and 1. Modify in',path)
		raise

	#Check that radius is an integer
	if type(D['radius']) == int:
		radius = D['radius']
	else:
		print('Radius input must be an integer. Modify in',path)
		raise

	#Check genthresh float between 0 and 1
	if D['genthresh'] <= 1 and D['genthresh'] >= 0:
		genthresh = D['genthresh']
	else:
		print('General threshold input must be a number between 0 and 1. Modify in',path)
		raise

	#Check format of microns
	if type(D['microns']) == list:
		for i in D['microns']:
			#If it is a float or int
			if not (type(i) == float) | (type(i) == int):
				print('Micron inputs must be numeric values. Modify in',path)
				raise
		microns = D['microns']
	else:
		print('Micron input must be a list of numeric values. Modify in',path)
		raise

	#Check degree
	if type(D['deg']) != int:
		print('Degree input must be an integer. Modify in',path)
		raise

# #Check component order
# 		lv = []
# 		for key in ['entry_x','entry_y','entry_z']:
# 			s = p['comporder'][key].get()
# 			try:
# 				v = int(s)
# 				lv.append(v)
# 			except ValueError:
# 				messagebox.showerror('Error','Component order inputs must be numeric values')
# 				return
# 			pc['comporder'] = lv

# 		#Get fit dimensions
# 		s = p['fitdim']['entry'].get()
# 		if s == 'XZ Plane':
# 			pc['fitdim'] = ['x','z']
# 		elif s == 'XY Plane':
# 			pc['fitdim'] = ['x','y']
# 		elif s == 'YZ Plane':
# 			pc['fitdim'] = ['y','z']

# 		for key in pc.keys():
# 			print(key,pc[key])

# 		return(pc)	

if __name__=='__main__':

	config_path = argv
	config_data = open(config_path).read()
	params = json.loads(config_data)
	check_config(params,config_path)