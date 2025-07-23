########################################################

# Code and data for 
# ‘Temperature-driven scaling patterns emerge in diatom 
# gene expression across the global ocean'

########################################################

# R script that filters the meta-omics abundance data 
# per sample by the diatom taxonomy, the minimum
# abundance of 10 and the presence in metaG for metaT
# unigenes. The output is an abundance table with the 
# rows equal to unigenes and columns sampling stations.

########################################################

# Required packages:
library(tidyverse)

########################################################

# Read the original abundance data per sample for both
# metaG and metaT.
list_files_G <- list.files(path='metaG_micro/by_sampleID',full.names=TRUE)
list_files_T <- list.files(path='metaT_micro/by_sampleID',full.names=TRUE)

# Get the unigene MATOU ID for Bacillariophyta
ID_bacilla <- read.csv('ID_bacilla.csv',header=TRUE) %>% pull(geneID)

# Create the containers for the final abundance tables
metaG_tot <- matrix(nrow=0,ncol=3) %>% as.data.frame() %>%
				rename(geneID=V1,counts=V2,station=V3)
metaT_tot <- matrix(nrow=0,ncol=3) %>% as.data.frame() %>%
				rename(geneID=V1,counts=V2,station=V3)

# for each station in the data
for(i in 1:length(list_files_G)){

	# read the metaG data
	metaG <- read.csv(list_files_G[i],header=TRUE) %>% 
					# get the diatom unigenes
					filter(geneID %in% ID_bacilla) %>%
					# remove abundances lower than 10
					filter(readCount >= 10) %>%
					rename(countsG=readCount)
	# read the metaT data
	metaT <- read.csv(list_files_T[i],header=TRUE) %>% 
					filter(geneID %in% ID_bacilla) %>%
					filter(readCount >= 10) %>%
					rename(countsT=readCount)
	# get the station ID
	stat_i <- unlist(strsplit(unlist(strsplit(list_files_G[i],split="/"))[3],
						split="SUR"))[1] %>% as.numeric()
	print(stat_i)	

	# Get the intersection of metaG and metaT unigenes 
	# and use this to filter metaT
	metaT <- inner_join(metaG,metaT,by='geneID') %>%
				select(c(geneID,countsT))

	# Add the unigenes of this station to the container 
	metaG_tot <- rbind(metaG_tot, 
						metaG %>% rename(counts=countsG) %>% 
									mutate(station=stat_i))
	metaT_tot <- rbind(metaT_tot, 
						metaT %>% rename(counts=countsT) %>% 
									mutate(station=stat_i))

}

# Transform the final data tables to unigene x station
# format
metaG_tot <- metaG_tot %>% pivot_wider(names_from=station,values_from=counts)
metaT_tot <- metaT_tot %>% pivot_wider(names_from=station,values_from=counts)

# Reorder the columns/stations
cols <- names(metaT_tot)[-1] %>% as.numeric() %>% sort() %>% as.character()
cols <- c('geneID',cols)

# Write the final tables to CSV files
metaG_tot[,cols] %>% 
			write.csv('metaG_micro_bacilla.csv',row.names=FALSE)
metaT_tot[,cols] %>% 
			write.csv('metaT_micro_bacilla.csv',row.names=FALSE)
