from cranium import landmark

#Create a optimization object
opt = landmark.anumSelect(dfs)

tstep = np.pi/4

#Initiate parameter sweep
opt.param_sweep(tstep,amn=2,amx=50,step=1,percbins=[50],rnull=15)

#Plot raw data
opt.plot_rawdata()

poly_degree = 4

#Test polynomial fit
opt.plot_fitted(poly_degree)

best_guess = 30

#Find the optimum value of anum
opt.find_optimum_anum(poly_degree,best_guess)
