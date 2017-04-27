import cranium
import multiprocessing as mp
import time
import os
import sys
from pushbullet import Pushbullet
#import plotly

#plotly.tools.set_credentials_file(username='msschwartz21', api_key='OrM0hMDBvseeT6SCjxNb')


if __name__=='__main__':

	dirs = ["C:\\Users\\zfishlab\\Desktop\\zrf1wt13umyot21um\\wt\\AT2",
	"C:\\Users\\zfishlab\\Desktop\\zrf1wt13umyot21um\\wt\\ZRF1",
	"C:\\Users\\zfishlab\\Desktop\\zrf1wt13umyot21um\\yot\\AT",
	"C:\\Users\\zfishlab\\Desktop\\zrf1wt13umyot21um\\yot\\ZRF1"]

	for d in dirs:
		print(d)
		files = os.listdir(d)
		
		#Check that all files are hdf5
		hfiles = []
		for file in files:
			if 'h5' in file and 'Probabilities' in file:
				hfiles.append(os.path.join(d,file))

		pool = mp.Pool()
		pool.map(cranium.process_sample,hfiles)
		pool.close()
		pool.join()
