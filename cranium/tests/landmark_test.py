from cranium import landmark
import numpy as np
import pandas as pd

# TODO add psi data to sample data directory

anum = 30
tstep = np.pi/4

#Create a landmark object
lm = landmark.landmarks(percbins=[50],rnull=np.nan)
lm.calc_bins(dfs,anum,tstep)

#Calculate landmarks for each sample and append to a single dataframe
outlm = pd.DataFrame()
for k in dfs.keys():
    outlm = lm.calc_perc(dfs[k],k,'stype',outlm)
