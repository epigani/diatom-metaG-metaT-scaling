########################################################

# Code and data for 
# ‘Temperature-driven scaling patterns emerge in diatom 
# gene expression across the global ocean'

########################################################

# Python script that calculates the minimum abundance 
# value μ by minimising the Kolmogorov-Smirnov (KS) 
# statistic of the fit of the theoretical distribution 
# to the filtered data.

# Input parameter (line 176):
# metaGT = metaG -> perform the analysis on the 
#	metagenomics data
# metaGT = metaT -> perform the analysis on the 
#	metatranscriptomics data

########################################################

# Required modules:
import numpy as np
import scipy as sp
import pandas as pd

########################################################

# Helper function with the theoretical distribution (PDF)
# This is the discretised version of the original 
# continuous distribution.
# - n: vector of unigene abundances
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def pareto3(n,a,logk,mu):

	# use the difference between the CDF 
	# at (n+1) with the one at (n)
	# to discretise the distribution
	return( pow((pow(10,logk) + n - mu)/pow(10,logk),-a) - pow((pow(10,logk) + n + 1.0 - mu)/pow(10,logk),-a) )

########################################################

# Helper function that returns the loglikelihood 
# function of the theoretical distribution using the 
# abundance values
# - fit_args: vector with the parameters of the 
#				theoretical distribution: [a,logk]
# - n_array: vector of unigene abundances
def loglik_par(fit_args,n_array):

	# get the individual parameter values
	a = fit_args[0]
	logk = fit_args[1]
	mu = np.min(n_array)

	# definition of the loglikelihood
	# minus sign for minimisation
	return(-np.sum(np.log(pareto3(n_array,a,logk,mu))))

########################################################

# Helper function that minimises the loglikelihood 
# function to find the optimal parameter values
# - n_array: vector of unigene abundances
# - fit_init: initial seed parameter values
def minimise_LL_par(n_array,fit_init):

	# the α value cannot be negative (0,None) 
	res = sp.optimize.minimize(fun=loglik_par,x0=fit_init,args=n_array,bounds=((0,None),(None,None)),method="Nelder-Mead")
	# return the fit results
	return(res)

########################################################

# Helper function with the theoretical cumulative 
# distribution function (CDF)
# - n: vector of unigene abundances
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def cumul(n,a,logk,mu):

	return(1.0 - pow((pow(10,logk)+n+1.0-mu)/pow(10,logk),-a))

########################################################

# Helper function that calculates the Kolmogorov-Smirnov
# statistic of the fit of the theoretical distribution
# to the filtered data
# - mu: the minimum abundance value μ
# - n_array: vector of unigene abundances
# - fit_init: initial seed parameter values
def KS(mu,n_array,fit_init):

	# filter the unigene abundances (minimum abundance)
	n_array = list(filter(lambda num: num >= mu, n_array))

	# fit the theoretical distribution to the filtered 
	# data and extract the individual parameters
	res = minimise_LL_par(n_array,fit_init).x
	a = res[0]
	logk = res[1]

	# construct the empirical CDF
	counts = np.asarray(pd.Series(n_array).value_counts().sort_index(ascending=True))
	d = {'x':np.unique(n_array),'S':np.cumsum(counts)/np.sum(counts)}

	# create a dataframe
	df_cumul = pd.DataFrame(data=d)
	# calculate the theoretical CDF
	df_cumul['P'] = cumul(df_cumul['x'].values,a,logk,mu)
	# return the KS statistic = maximum difference 
	# between the theoretical and empirical CDF
	return(np.max(abs(df_cumul['S'].values-df_cumul['P'].values)))

########################################################

# Helper function that finds the minimum abundance value
# that minimises the KS statistic of the resulting fit
# - n_array: vector of unigene abundances
# - fit_init: initial seed parameter values
# - mu_upper: the maximum μ value to test for 
#				minimising the KS statistic
def minimise_KS(n_array,fit_init,mu_upper):

	# create a vector of μ values to test
	th = np.arange(10,mu_upper)
	# create a vector to store the KS 
	# statistic of the resulting fit for
	# each tested μ value
	KS_th = np.zeros(len(th))
	# calculate the KS statistic for each
	# μ value
	for i in range(0,len(th)):
		KS_th[i] = KS(th[i],n_array,fit_init)

	# return the μ value that minimises the
	# KS statistic
	return(th[np.argmin(KS_th)])

########################################################

# Helper function that calculates the minimum abundance
# value μ by discarding the abundaces below the abundance 
# value for which the CDF is equal to the threshold value 
# (cum_thresh)
# - n_array: vector of unigene abundances
# - cum_thresh: the threshold value of the CDF
def get_cumthresh(n_array,cum_thresh):

	# the total abundance
	N_tot = sum(n_array)
	# construct the empirical CDF
	n_array = n_array/N_tot
	n_array = np.asarray(sorted(n_array))
	n_array_cumsum = np.asarray(np.cumsum(n_array))
	n_array = n_array*N_tot
	# find the abundance at which the CDF
	# surpasses the threshold
	thresh_val = n_array[(len(n_array)-len(list(filter(lambda num: num > cum_thresh, n_array_cumsum))))]

	return(thresh_val)

########################################################

# Main execution of the calculations
def main():

	# the input parameter
	# 'metaG' = do for metagenomics data
	# 'metaT' = do for metatranscriptomics data
	metaGT = "metaG"

	# read the unigene abundance table
	# set NA-values to 0
	meta = pd.read_csv('../data/'+metaGT+'_micro_bacilla.csv').fillna(0)
	# read the maximum μ values to test
	limits = pd.read_csv("../data/limit_KS_threshold_"+metaGT+"_micro_bacilla.csv")

	# get a list of Tara station IDs (string)
	stats = meta.columns
	stats = stats[1:len(stats)]

	# do things seperately for the first station to
	# correctly write to the CSV file
	print("station "+str(1)+" of "+str(len(stats)))
		
	# get the vector of unigene abundances and remove
	# absent unigenes 	
	n_array = meta[stats[0]].values
	n_array = list(filter(lambda num: num != 0, n_array))

	# get the μ value by using the threshold of the CDF
	cum_th = get_cumthresh(n_array,0.01)

	# get the maximum μ value to test for this station
	mu_upper = limits[limits["station"]==int(stats[0])]["mu_upper"][0]

	# default initial seed values for the parameters
	# (the final result is not sensitive to this choice)
	fit_init = [1.0,np.log10(cum_th)]

	# initialise the output CSV file with its 3 columns
	# by writing the output for the first station	
	pd.DataFrame(data={"station":[int(stats[0])],"mu_cum":[int(cum_th)],"mu_KS":[int(minimise_KS(n_array,fit_init,mu_upper))]}).to_csv('summary_KS_threshold_'+metaGT+'_micro_bacilla.csv',index=False)

	# do the same for all other stations and append
	# the output to the output CSV file
	for i in range(1,len(stats)):

		print("station "+str(i+1)+" of "+str(len(stats)))
		
		n_array = meta[stats[i]].values
		n_array = list(filter(lambda num: num != 0, n_array))

		cum_th = get_cumthresh(n_array,0.01)
		mu_upper = limits[limits["station"]==int(stats[i])]["mu_upper"][i]

		fit_init = [1.0,2.0]
		
		pd.DataFrame(data={"station":[int(stats[i])],"mu_cum":[int(cum_th)],"mu_KS":[int(minimise_KS(n_array,fit_init,mu_upper))]}).to_csv('summary_KS_threshold_'+metaGT+'_micro_bacilla.csv',index=False,header=False,mode='a')

main()