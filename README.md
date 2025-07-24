# diatom-metaG-metaT-scaling
Code and data for ‘Temperature-driven scaling patterns emerge in diatom gene expression across the global ocean'

## data

### filter

This directory contains the scripts and tables that are used in processing and filtering the original [MATOU](https://www.genoscope.cns.fr/tara/) unigene abundance tables.  
The downloaded original abundance table is first filtered by size class (20-180 μm) and sampling depth (SURFACE). This very large data file is cut into smaller pieces so that they can be loaded in memory. These files are characterised by:  

*columns:*

* 'geneID': MATOU unigene ID;
* 'sampleName': MATOU sample ID;
* 'readCount': abundance of the unigene in the sample.

*rows:*  
individual data entries (unigene-sample)  
  
The Python script **preprocessing.py** takes these raw abundance files and creates abundance files per station. These files contain:  

*columns:*

* 'geneID': MATOU unigene ID;
* 'readCount': abundance of the unigene.

*rows:*  
individual unigenes  
  
The R script **get_taxon.R** takes the MATOU taxonomy table. This taxonomy table lists the taxonomy of each unigene in the MATOU data:  

*columns:*

* 'geneID': MATOU unigene ID;
* 'taxID': MATOU taxon ID;
* 'taxName': name of the lowest taxonomic level;
* 'taxRank': the lowest taxonomic level;
* 'taxLineage': the full taxonomy up to the lowest level.

*rows:*  
individual unigenes (76,133,096)

The script extracts all the unigene IDs that belong to the Bacillariophyta taxon. These unigene IDs are listed in the **ID_bacilla.csv** table under the 'geneID' column.  
  
Finally, the **filter_counts.R** script filters the original unigene abundance data for metaG and metaT using the Bacillariophyta unigene IDs, discarding all abundances lower than 10 and taking the intersection with metaG for the metaT data.  
  
The results are the **metaG_micro_bacilla.csv** and **metaT_micro_bacilla.csv** tables. These filtered MATOU unigene abundance tables can be found in the **metaGT_micro_bacilla.zip** file. Before executing the fitting scipts, first extract the abundance tables from the zip file into the **data/filter** directory.  
  
**metaG_micro_bacilla.csv**: abundance table for metagenomics unigenes filtered on the 20-180 μm (micro) size fraction and the diatom (Bacillariophyta) class.  
**metaT_micro_bacilla.csv**: abundance table for metatranscriptomics unigenes filtered on the 20-180 μm (micro) size fraction and the diatom (Bacillariophyta) class.  
  
*columns*:

* 'geneID': MATOU unigene ID; 
* other columns: Tara Oceans station ID.

*rows*:  
individual unigenes  

*values*:

* abundance of unigene X in station Y;
* fill value: NA.

### environment

The **woa_pisces_tarassd.tsv** table contains environmental information about a subset of Tara Oceans surface stations. In the analysis, only the 'temp_woa' column is used for the average Sea Surface Temperature at the sampling locations.

## fit

### KS threshold parameters

To determine the minimum abundance value μ, unigenes are progressively removed starting from the lowest abundances. The value μ for which the Kolmogorov-Smirnov (KS) statistic of the fit to the remaining abundance data is minimised is taken to be the minimum abundance value for the following fitting procedure. Here, the upper μ value to be tested in this way is determined by visually inspecting the cumulative distribution function (CDF) to find a balance between computing time and accuracy. More automated methods were tested but do not significantly change the final results while making the computation impractical.

**limit_KS_threshold_metaG_micro_bacilla.csv**: Helper table for the maximum μ values to be tested for each station for the metagenomics filtered data set.  
**limit_KS_threshold_metaT_micro_bacilla.csv**: Helper table for the maximum μ values to be tested for each station for the metatranscriptomics filtered data set.  
  
*columns*:

* 'station': Tara Oceans station ID; 
* 'mu_upper': maximum μ value to be tested.

*rows*:  
individual stations

### minimum abundance values

Results of the algorithm to minimise the KS statistic by progressively increasing the minimum abundance value described above. This algorithm is implemented in the *calc_KS_threshold.py* script.  
  
**summary_KS_threshold_metaG_micro_bacilla.csv**: minimum abundance value used in the fitting procedure for each station for the metagenomics filtered data set.  
**summary_KS_threshold_metaT_micro_bacilla.csv**: minimum abundance value used in the fitting procedure for each station for the metatranscriptomics filtered data set.  

*columns*:

* 'station': Tara Oceans station ID; 
* 'mu_cum': minimum abundance value μ determined by taking the abundance where the CDF is 0.01 (not used in further analyses);
* 'mu_KS': minimum abundance value μ determined by minimising the KS statistic as described above.

*rows*:  
individual stations

### fit results

The following fitting procedure is applied to each station. The filtered abundances are further filtered by removing all unigenes with abundance lower than the minimal abundance μ found before. Next, the theoretical curve is fitted to the original abundances. For the goodness-of-fit test, replicate abundances are drawn from the fitted theoretical distribution. The theoretical distribution is then fitted to these replicate samples and the statistics are recalculated. The final statistics are calculated from the combined samples. The Python script *fit_pareto.py* performs this algorithm.  
  
**summary_fit_metaG_micro_bacilla.csv**: summary of the fitting procedure and goodness-of-fit test for each station for the filtered metagenomics data set.  
**summary_fit_metaT_micro_bacilla.csv**: summary of the fitting procedure and goodness-of-fit test for each station for the filtered metatranscriptomics data set.  

*columns*:

* 'station': Tara Oceans station ID; 
* 'mu_KS': minimum abundance value μ determined by minimising the KS statistic as described above;  
* 'N_10': total abundance of unigene reads for the initially filtered data;
* 'S_10': number of unique unigenes for the initially filtered data;  
* 'N_KS': total abundance of unigene reads filtered by the KS minimum threshold;
* 'S_KS': number of unique unigenes filtered by the KS minimum threshold;
* 'alpha': the α parameter of the theoretical distribution fitted to the original filtered data;
* 'alpha_low': the 2.5-percent quantile of the α parameter of the theoretical distribution fitted to the replicates;
* 'alpha_high': the 97.5-percent quantile of the α parameter of the theoretical distribution fitted to the replicates;
* 'logk': the 10-based logarithm of the k parameter of the theoretical distribution fitted to the original filtered data;
* 'logk_low': the 2.5-percent quantile of the 10-based logarithm of the k parameter of the theoretical distribution fitted to the replicates;
* 'logk_high': the 97.5-percent quantile of the 10-based logarithm of the k parameter of the theoretical distribution fitted to the replicates;
* 'logN_mean': the mean of the 10-based logarithm of the total abundance of unigene reads in the replicates;
* 'logN_low': the 2.5-percent quantile of the 10-based logarithm of the total abundance of unigene reads in the replicates;
* 'logN_high': the 97.5-percent quantile of the 10-based logarithm of the total abundance of unigene reads in the replicates;
* 'KS': the Kolmogorov-Smirnov statistic of the fit of the theoretical distribuution to the original filtered data;
* 'p_KS': the fraction of the replicates with a KS statistic higher than the one of the fit to the original filtered data;
* 'AD': the Anderson-Darling statistic of the fit of the theoretical distribuution to the original filtered data;
* 'p_AD': the fraction of the replicates with an AD statistic higher than the one of the fit to the original filtered data;


*rows*:  
individual stations

### scripts

**fit/calc_KS_threshold.py**: Python script that calculates the minimum abundance value μ by minimising the Kolmogorov-Smirnov (KS) statistic of the fit of the theoretical distribution to the filtered data.  
The only input parameter is whether to perform the analysis on the metagenomics (metaG) or metatranscriptomics (metaT) data and can be changed directly in the main function of the script.  
The script reads the required abundance table and the table with the KS threshold values.  
The output of the script is the table with the minimum abundance values.  
The programme takes around 2 hours to run.  
  
**fit/fit_pareto.py**: Python script that performs the fit and goodness-of-fit test of the theoretical distribution to the filtered abundance tables.  
In the main function of the script, two parameters can be adjusted. The first input parameter is whether to perform the analysis on the metagenomics (metaG) or metatranscriptomics (metaT) data. The second parameter is the number of replicates to use for the goodness-of-fit test. If Nrep is the number of replicates, the minimum p-value for the test is 1/Nrep.  
The output of the script is a table with the summary of the fit of the theoretical distribution and the goodness-of-fit test for all stations.  
The programme takes a few hours to run for 1,000 replicates.