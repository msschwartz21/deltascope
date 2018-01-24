import os
import re
import pandas as pd

def write_header(f):
	'''
	Writes header for PSI file with columns Id,x,y,z,ac,r,theta

	:param file f: file object created by 'open(filename,'w')`
	'''

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
	'''
	Writes data in PSI format to file after writing header using :py:func:`write_header`. Closes file at the conclusion of writing data.

	:param str filepath: Complete filepath to output file
	:param pd.DataFrame df: dataframe containing columns x,y,z,ac,r,theta
	'''

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
	'''
	Reads psi file at the given filepath and returns data in a pandas DataFrame

	:param str filepath: Complete filepath to file
	:returns: pd.Dataframe containing data
	'''

	df = pd.read_csv(filepath,
		sep=' ',
		header=19)

	if len(df.columns) <= 4:
		df.columns = ['i','x','y','z']
	elif len(df.columns) >= 6:
		df.columns=['i','x','y','z','ac','theta','r']

	return(df)

def read_psi_to_dict(directory,dtype):
	'''
	Read psis from directory into dictionary of dfs with filtering based on dtype

	:param str directory: Directory to get psis from
	:param str dtype: Usually 'AT' or 'ZRF1'
	:returns: Dictionary of pd.DataFrame
	'''

	dfs = {}
	for f in os.listdir(directory):
		if dtype in f:
			df = read_psi(os.path.join(directory,f))
			num = re.findall(r'\d+',f.split('.')[0])[0]
			print(num)
			dfs[num] = df

	return(dfs)
