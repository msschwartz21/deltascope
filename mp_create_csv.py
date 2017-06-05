import cranium
import multiprocessing as mp
import time
import os
import sys
from pushbullet import Pushbullet
from functools import partial
from random import randint
#import plotly

#plotly.tools.set_credentials_file(username='msschwartz21', api_key='OrM0hMDBvseeT6SCjxNb')


if __name__=='__main__':

	root = 'C:\\Users\\zfishlab\\Desktop\\zrf1wt13umyot21um'
	stypes = ['yot'] #['wt','yot']
	chs = ['AT2\\Prob','ZRF1']
	prefixes = ['AT','ZRF1']

	#Create outdir
	outdir = os.path.join(root,'Output'+str(randint(10,99)))
	os.mkdir(outdir)
	print('outdir',outdir)

	for stype in stypes:
		d = os.path.join(root,stype,chs[0])
		print(stype,'starting')

		files = os.listdir(d)
		#Get sample numbers
		nums = []
		for file in files:
			if 'h5' in file and 'Probabilities' in file:
				nums.append(file.split('_')[1])

		print(nums)

		processfxn = partial(cranium.process_sample,
			root=os.path.join(root,stype),
			outdir=outdir,
			name=stype,
			chs=chs,
			prefixes=prefixes,
			threshold=0.5,
			scale=[1000,100,1],
			deg=2,
			primary_key='at')

		n = len(nums)
		for i in range(0,n,5):
			if i+5 > n:
				L = nums[i:n]
			else:
				L = nums[i:i+5]

			pool = mp.Pool()
			pool.map(processfxn,L)
			pool.close()
			pool.join()

		print(stype, 'complete')


	# for d in dirs:
	# 	print(d)
	# 	files = os.listdir(d)
		
	# 	#Check that all files are hdf5
	# 	hfiles = []
	# 	for file in files:
	# 		if 'h5' in file and 'Probabilities' in file:
	# 			hfiles.append(os.path.join(d,file))

	# 	pool = mp.Pool()
	# 	pool.map(cranium.process_sample,hfiles)
	# 	pool.close()
	# 	pool.join()

	# 	print(d,'complete')
