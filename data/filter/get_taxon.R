########################################################

# Code and data for 
# ‘Temperature-driven scaling patterns emerge in diatom 
# gene expression across the global ocean'

########################################################

# R script that extracts the MATOU gene IDs of the 
# unigenes belonging to the class Bacillariophyta. 

########################################################

# Required packages
library(tidyverse)
library(data.table)

########################################################

# Read the MATOU taxonomy table that contains the 
# taxonomy of each MATOU unigene
# This table contains 76,133,096 unique unigenes.
dat <- read.table('MATOU-v2.taxonomy.tsv',header=TRUE,sep='\t')

# the diatom taxon
taxon <- "Bacillariophyta"
taxon_ID <- "bacilla"

# Divide the full taxonomy table in pieces of 100,000
# unigenes to fit them into memory.
# Do the first 100,000 seperately to create the output
# file with diatom gene IDs.

# Deflate the taxLineage column, which contains the full
# taxonomy of the unigene
# According to this seperation, the Bacillariophyta level
# will be a 'subkingdom'.
test <- separate_wider_delim(dat[((0)*1e5+1):((1)*1e5),], cols = taxLineage, delim = ";",
					names=c('root','group','domain','clade','kingdom','subkingdom','group_2',
								'group_3','supergroup','superphylum','phylum','subphylum_1','subphylum_2',
								'subphylum_3','subphylum_4','superclass','class','order','family',
								'genus','species'),
					too_few="align_start",
					too_many='drop')

# Get the MATOU gene IDs for all unigenes belonging to
# Bacillariophyta. 
vec_ID <- test %>% dplyr::filter(subkingdom==taxon | taxName==taxon) %>% pull(geneID)

# Create the output file with gene IDs and write the IDs
# from the first 100,000 lines.
write.table(data.frame(geneID=vec_ID),paste0('ID_',taxon_ID,'.csv'),sep=',',col.names=TRUE,append=FALSE,row.names=FALSE)

# Now do the same for all other pieces of 100,000 rows.
for(i in 1:760){

	print(i)

	test <- separate_wider_delim(dat[((i)*1e5+1):((i+1)*1e5),], cols = taxLineage, delim = ";",
					names=c('root','group','domain','clade','kingdom','subkingdom','group_2',
								'group_3','supergroup','superphylum','phylum','subphylum_1','subphylum_2',
								'subphylum_3','subphylum_4','superclass','class','order','family',
								'genus','species'),
					too_few="align_start",
					too_many='drop')

	vec_ID <- test %>% dplyr::filter(subkingdom==taxon | taxName==taxon) %>% pull(geneID)

	# add these diatom unigenes to the output file
	write.table(data.frame(geneID=vec_ID),paste0('ID_',taxon_ID,'.csv'),sep=',',col.names=FALSE,append=TRUE,row.names=FALSE)
}

# Finally do the same for the last 100,000 unigenes
test <- separate_wider_delim(dat[(76100001):(76133096),], cols = taxLineage, delim = ";",
					names=c('root','group','domain','clade','kingdom','subkingdom','group_2',
								'group_3','supergroup','superphylum','phylum','subphylum_1','subphylum_2',
								'subphylum_3','subphylum_4','superclass','class','order','family',
								'genus','species'),
					too_few="align_start",
					too_many='drop')

vec_ID <- test %>% dplyr::filter(subkingdom==taxon | taxName==taxon) %>% pull(geneID)

write.table(data.frame(geneID=vec_ID),paste0('ID_',taxon_ID,'.csv'),sep=',',col.names=FALSE,append=TRUE,row.names=FALSE)