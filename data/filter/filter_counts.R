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
	stat_i <- unlist(strsplit(unlist(strsplit(list_files_T[i],split="/"))[3],
						split="SUR"))[1] %>% as.numeric()
	print(stat_i)
	

	# Get the intersection of metaG and metaT unigenes 
	# and use this for filtering metaT
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
	
	if(F){
	stat_all <- c( 7,11,18,22,23,25,26,30,36,38,39,41,46,51,52,64,65,66,67,68,70,72,76,78,80,81,82,84,  
					85,92,98,100,102,109,110,111,122,123,124,125,128,131,132,135,136,137,138,142,143,144, 
					145,146,147,148,149,150,151,152,155,158,168,173,175,178,180,188,189,191,193,194,196, 
					201,205,206,208,209,210)

	stat_not <- stat_all[which(!(stat_all %in% (metaT_tot %>% pull(station) %>% unique()) ))]
	metaT_tot <- rbind(metaT_tot,
						data.frame(geneID=rep(metaT_tot$geneID[1],length(stat_not)),
									counts=rep(NA,length(stat_not)),
									station=stat_not))
	stat_not <- stat_all[which(!(stat_all %in% (metaG_tot %>% pull(station) %>% unique()) ))]
	metaG_tot <- rbind(metaG_tot,
						data.frame(geneID=rep(metaG_tot$geneID[1],length(stat_not)),
									counts=rep(NA,length(stat_not)),
									station=stat_not))

	}

	metaG_tot <- metaG_tot %>% pivot_wider(names_from=station,values_from=counts)
	metaT_tot <- metaT_tot %>% pivot_wider(names_from=station,values_from=counts)

	cols <- names(metaT_tot)[-1] %>% as.numeric() %>% sort() %>% as.character()
	cols <- c('geneID',cols)

	metaG_tot[,cols] %>% 
				write.csv('metaG_micro_bacilla.csv',row.names=FALSE)
	metaT_tot[,cols] %>% 
				write.csv('metaT_micro_bacilla.csv',row.names=FALSE)
