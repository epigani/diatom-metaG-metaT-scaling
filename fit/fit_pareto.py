########################################################

# Code and data for 
# ‘Temperature-driven scaling patterns emerge in diatom 
# gene expression across the global ocean'

########################################################

# Python script that fits the theoretical distribution
# to the filtered Tara Oceans unigene abundances and 
# performs the goodness-of-fit tests.

# Input parameters (line 313):
# metaGT = metaG -> perform the analysis on the 
#	metagenomics data
# metaGT = metaT -> perform the analysis on the 
#	metatranscriptomics data
# Nrep = the number of replicates for the goodness-of-fit
# 			tests

########################################################

# Required modules:
import scipy as sp
import numpy as np
import pandas as pd
import random

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

# Helper function with the inverse of the theoretical CDF.
# This expression can be analytically derived.
# - r: a number drawn from a uniform distribution 
#		between 0 and 1
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def cumul_inv(r,a,logk,mu):

	return( pow(10,logk)*pow(1.0-r,-1.0/a) - pow(10,logk) + mu - 1.0 )

########################################################

# Helper function that samples a random number from the 
# theoretical distribution by using the analytical 
# expression of the inverse of the theoretical CDF
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def sample_pareto_inv(a,logk,mu):

	# draw from the uniform distribution
	r = random.uniform(0,1)
	# always round up for the correct minimum abundance
	# value
	return( np.ceil(cumul_inv(r,a,logk,mu)) )

########################################################

# Helper function that generates a full sample from the
# fitted theoretical distribution
# - N: the number of unigene abundances to draw
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def create_sub(N,a,logk,mu):

	# vector for the sampled abundances
	out = np.zeros(N)
	# fill the vector
	for i in range(0,N):
		out[i] = sample_pareto_inv(a,logk,mu)

	return(out)

########################################################

# Helper function that calculates the Kolmogorov-Smirnov
# statistic of the fit of the theoretical distribution
# to the abundance vector
# - n_array: vector of unigene abundances
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def KS(n_array,a,logk,mu):

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

# Helper function that calculates the Anderson-Darling
# statistic of the fit of the theoretical distribution
# to the abundance vector
# - n_array: vector of unigene abundances
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def AD(n_array,a,logk,mu):

	# initialise the S sum
	S_stat = 0
	# sort the abundances
	n_array = np.asarray(sorted(n_array))
	# calculate the S sum
	for i in range(1,len(n_array)+1):
		S_stat += (2*i-1)*(np.log(cumul(n_array[i-1],a,logk,mu)) + np.log(1.0-cumul(n_array[len(n_array)-i],a,logk,mu)))
	# return the AD statistic by correcting the S sum
	return(np.sqrt(-len(n_array)-S_stat/len(n_array)))

########################################################

# Helper function that calculates the empirical evenness
# of the abundance vector
# This measure is not used in further analyses.
# - n_array: vector of unigene abundances
def evenness(n_array):

	# calculate relative abundances
	p_array = n_array/np.sum(n_array)
	# return the Shannon index
	return(-np.sum(p_array*np.log(p_array))/np.log(len(n_array)))

########################################################

# Helper function that calculates the statistics
# - n_array: vector of unigene abundances
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def calc_statistics(n_array,a,logk,mu):

	# the goodness-of-fit statistsics
	# Kolmogorov-Smirnov
	KS_ori = KS(n_array,a,logk,mu)
	# Anderson-Darling
	AD_ori = AD(n_array,a,logk,mu)

	# calculate the evenness
	even = evenness(n_array)

	# return some interesting statistics as a vector
	return([np.sum(n_array),a,logk,KS_ori,AD_ori,even,np.median(n_array),np.mean(n_array)])

########################################################

# Helper function that generates the replicate sample and
# fits the theoretical distribution to the replicate
# - S: the number of unigene abundances in the sample
# - a: α parameter of the theoretical distribution
# - logk: 10-based logarithm of the k parameter of the 
#			theoretical distribution
# - mu: the minimum abundance value μ
def replicate_run(S,a,logk,mu):

	# create the replicate abundances
	sub_array = create_sub(S,a,logk,mu)
	# calculate the replicate total abundance
	N_sub = np.sum(sub_array)
	# fit the theoretical distribution to the 
	# replicate sample
	res = minimise_LL_par(sub_array,[a,logk]).x
	# get the individual parameters
	a_sub = res[0]
	logk_sub = res[1]

	# return some interesting statistics as a vector
	return(calc_statistics(sub_array,a_sub,logk_sub,mu))

########################################################

# Helper function that performs the fit and 
# goodness-of-fit tests for a station:
# 1) fit the theoretical distribution to the original 
#		filtered abundances
# 2) replicate the sample Nrep times and do the same as in 1)
# 3) calculate the summary statistics
# - stat_ID: the Tara Oceans station ID
# - n_array: vector of original filtered unigene abundances
# - mu: the minimum abundance value μ 
# - Nrep: number of replicates to perform
def do_station(stat_ID,n_array,mu,Nrep):

	# store the original abundance vector
	full_array = n_array
	# discard all abundances lower than the threshold
	n_array = n_array[n_array >= mu]
	# get the number of unique unigenes remaining in
	# the sample
	S = len(n_array)

	# default initial seed values for the parameters
	# (the final result is not sensitive to this choice)
	fit_init = [1.0,2.0]

	# container for the statistics of the original and
	# replicate calculations
	replicates = np.zeros((Nrep,8))

	# fit the theoretical distribution to the original
	# filtered abundance data
	res = minimise_LL_par(n_array,fit_init).x

	# the first row is the fit to the original data
	replicates[0,:] = calc_statistics(n_array,res[0],res[1],mu)

	# the remaining rows are the fits to the 
	# replicate data runs
	for i in range(1,Nrep):
		replicates[i,:] = replicate_run(S,res[0],res[1],mu)

	# return the summary statistics as described
	# in the README file
	return([stat_ID,mu,np.sum(full_array),len(full_array[full_array > 0]),
				np.sum(n_array),S,replicates[0,1],
				np.quantile(replicates[:,1],0.025),
				np.quantile(replicates[:,1],0.975),
				replicates[0,2],np.quantile(replicates[:,2],0.025),
				np.quantile(replicates[:,2],0.975),
				np.mean(np.log10(replicates[:,0])),
				np.quantile(np.log10(replicates[:,0]),0.025),
				np.quantile(np.log10(replicates[:,0]),0.975),
				replicates[0,3],
				len(replicates[:,3][replicates[:,3] >= replicates[0,3]])/Nrep,
				replicates[0,4],
				len(replicates[:,4][replicates[:,4] >= replicates[0,4]])/Nrep])

########################################################

# Main execution of the calculations
def main():

	# the input parameters:
	# 'metaG' = do for metagenomics data
	# 'metaT' = do for metatranscriptomics data
	metaGT = "metaG"
	# number of replicates for the goodness-of-fit test
	Nrep = 1000

	# read the unigene abundance table
	# set NA-values to 0
	meta = pd.read_csv('../data/'+metaGT+'_micro_bacilla.csv').fillna(0)
	# read the minimum abundance values μ
	mu_values = pd.read_csv("../data/summary_KS_threshold_"+metaGT+"_micro_bacilla.csv")

	# get a list of Tara station IDs (string)
	stats = meta.columns
	stats = stats[1:len(stats)]

	# container for the summary statistics of all the
	# stations
	out = np.zeros((len(stats),19))
	
	# calculate the summary statistics for all stations
	for i in range(0,len(stats)):

		print("station "+str(i+1)+" of "+str(len(stats)))
		
		# get the vector of unigene abundances
		n_array = np.asarray(meta[stats[i]].values)

		# get the minimum abundance value μ
		mu = mu_values[mu_values["station"]==int(stats[i])]["mu_KS"][i]

		# store the summary statistics as a row
		out[i,:] = do_station(int(stats[i]),n_array,mu,Nrep)

	# create the output dataframe with the column names
	# as described in the README file
	df_out = pd.DataFrame(data=out,
		columns=["station","mu_KS","N_10","S_10","N_KS","S_KS",
					"alpha","alpha_low","alpha_high",
					"logk","logk_low","logk_high",
					"logN_mean","logN_low","logN_high",
					"KS","p_KS","AD","p_AD"])

	# change some columns to integer values
	df_out = df_out.astype({"station":np.int64,"mu_KS":np.int64,
				"N_10":np.int64,"S_10":np.int64,"N_KS":np.int64,
				"S_KS":np.int64})
	
	# write the dataframe to a CSV file
	df_out.to_csv("summary_fit_"+metaGT+"_micro_bacilla.csv",index=False)

main()