########################################################

# Code and data for 
# ‘Temperature-driven scaling patterns emerge in diatom 
# gene expression across the global ocean'

########################################################

# R script that calculates the mean Sea Surface 
# Temperatures at the Tara sampling stations by bilinear
# interpolation from World Ocean Atlas (WOA) data.

########################################################

# Required packages:
library(tidyverse)
library(ncdf4)
library(data.table)

########################################################

# Helper function that transforms WGS84 latitudes to 
# new coordinates that make the transformed grid 
# equidistant
# - lat: WGS84 latitude 
mulat <- function(lat){

	return(log(abs(1/cos(lat)+tan(lat))))

}

########################################################

# Helper function that performs the bilinear
# interpolation for a sampling station
# - x: longitude (radians) of the sampling station
# - y: mu coordinate of the sampling station
# - lons: vector of longitude (radians) values of the 
#			data grid
# - lats: vector of mu-coordinate values of the data grid
# - envmat: 2D SST grid as a 2D matrix
bilin <- function(x,y,lons,lats,envmat){

	# get the grid longitude to the left of the sample
	i_lon <- max(which(lons < x))
	x_lim <- c(lons[i_lon],lons[i_lon+1])
	# get the grid mu coordinate below the sample
	i_lat <- max(which(lats < y))
	y_lim <- c(lats[i_lat],lats[i_lat+1])
	
	# o   o
	#   x
	# o   o

	# get the horizontal distances of the square
	# around the sampling location
	x_mat <- matrix(nrow=1,ncol=2)
	x_mat[1,1] <- x_lim[2] - x
	x_mat[1,2] <- x - x_lim[1]

	# get the vertical distances of the square
	# around the sampling location
	y_mat <- matrix(nrow=2,ncol=1)
	y_mat[1,1] <- y_lim[2] - y
	y_mat[2,1] <- y - y_lim[1]

	# get the grid values of the square around
	# the sampling location
	Q_mat <- matrix(nrow=2,ncol=2)
	# lower left
	Q_mat[1,1] <- envmat[i_lon,i_lat]
	# upper left
	Q_mat[1,2] <- envmat[i_lon,i_lat+1]
	# lower right
	Q_mat[2,1] <- envmat[i_lon+1,i_lat]
	# upper right
	Q_mat[2,2] <- envmat[i_lon+1,i_lat+1]

	# if some of the grid values of the square 
	# are missing
	if(is.na(sum(Q_mat))){

		# get the missing values
		i_NA <- which(is.na(Q_mat),arr.ind=TRUE)
		for(j in 1:nrow(i_NA)){
			# set the missing values equal to
			# the mean of the non-missing values
			Q_mat[i_NA[j,1],i_NA[j,2]] <- mean(Q_mat,na.rm=TRUE)
		}

		# if all grid values of the square are
		# missing
		if(is.na(sum(Q_mat))){

			# expand the grid:
			# o       o
			#   .   .
			#     x
			#   .   .
			# o       o
			# and do the same as before:
			x_lim <- c(lons[i_lon-1],lons[i_lon+2])
			y_lim <- c(lats[i_lat-1],lats[i_lat+2])
			x_mat[1,1] <- x_lim[2] - x
			x_mat[1,2] <- x - x_lim[1]
			y_mat[1,1] <- y_lim[2] - y
			y_mat[2,1] <- y - y_lim[1]

			Q_mat[1,1] <- envmat[i_lon-1,i_lat-1]
			Q_mat[1,2] <- envmat[i_lon-1,i_lat+2]
			Q_mat[2,1] <- envmat[i_lon+2,i_lat-1]
			Q_mat[2,2] <- envmat[i_lon+2,i_lat+2]

			# if some grid values of the expanded square
			# are missing
			i_NA <- which(is.na(Q_mat),arr.ind=TRUE)
			for(j in 1:nrow(i_NA)){
				# set the missing values equal to
				# the mean of the non-missing values
				Q_mat[i_NA[j,1],i_NA[j,2]] <- mean(Q_mat,na.rm=TRUE)
			}

			# if all grid values of the expanded square
			# are missing
			if(is.na(sum(Q_mat))){

				# expand the square even further
				# (up to 10 neighbours) and take all
				# values within the square 
				df_sub <- envmat[(i_lon-10):(i_lon+10),(i_lat-10):(i_lat+10)] %>% as.data.frame()
				names(df_sub) <- lats[(i_lat-10):(i_lat+10)]

				# now order the non-missing grid values in
				# this expanded square according to the 
				# distance to the sampling location of the 
				# grid points
				df_sub <- cbind(data.frame(lon=lons[(i_lon-10):(i_lon+10)]),df_sub) %>%
							pivot_longer(2:ncol(.),names_to="lat",values_to="SST") %>%
							mutate(lat=as.numeric(lat)) %>%
							filter(!is.na(SST)) %>%
							mutate(dist = sqrt((x-lon)^2 + (y-lat)^2)) %>%
							arrange(dist) 
				# extract the nearest grid point to the 
				# sampling location and take that grid 
				# value as the final value
				# if the data is the standard deviation,
				# then ignore zeros
				if(df_sub$SST[1] == 0){
					i_val <- mean(df_sub$SST[1:min(which(df_sub$SST != 0))])
				} else{
					i_val <- df_sub$SST[1]
				}				
				Q_mat[,] <- i_val
			}
		}
	}

	# return the bilinear interpolation
	return( as.numeric(x_mat %*% Q_mat %*% y_mat)/(x_lim[2]-x_lim[1])/(y_lim[2]-y_lim[1]) )

}

########################################################

# Main execution of the data extraction

# read the coordinates of the sampling stations
df_stat <- fread("TARA_reg_stations.tab",sep="\t",header=TRUE) %>% 
			select(all_of(c(2,5,6)))
names(df_stat) <- c("station","lat","lon")
# transform the coordinates to the equidistant coordinates
df_stat <- df_stat %>% mutate(lon=lon/180*pi,mu=mulat(lat/180*pi))

# read the SST data as downloaded from WOA
d <- nc_open("woa23_decav_t00_04.nc")

# extract both the mean SST as the 
# standard deviation (not used in further analyses)
# get the data at the second depth level (5 metres)
# to avoid inconsistencies at the real surface
temp_mn <- ncvar_get(d,"t_mn")[,,2]
temp_sd <- ncvar_get(d,"t_sd")[,,2] 

# get the longitude and latitude coordinates of the
# grid points and transform
lons <- d %>% ncvar_get("lon") 
lons <- lons/180*pi
lats <- d %>% ncvar_get("lat") 
mus <- mulat(lats/180*pi)

# create the final container for the interpolated 
# SST values at each sampling station
df_env <- matrix(ncol=5,nrow=nrow(df_stat)) %>% as.data.frame() %>%
			rename(station=V1,lon=V2,lat=V3,SST_mn=V4,SST_sd=V5)

# for each sampling station
for(i in 1:nrow(df_stat)){

	# interpolate the mean SST
	mn_i <- bilin(df_stat$lon[i],df_stat$lat[i],lons,lats,temp_mn)
	# interpolate the standard deviation of the SST
	sd_i <- bilin(df_stat$lon[i],df_stat$lat[i],lons,lats,temp_sd)
	# get the Tara Oceans station ID
	stat_i <- substr(df_stat$station[i],6,8) %>% as.numeric()

	# add the data to the container
	df_env[i,] <- c(stat_i,df_stat$lon[i]/pi*180,df_stat$lat[i],mn_i,sd_i)

}

# write the final data to a CSV file
df_env %>% write.csv("TARA_env_table.csv",row.names=FALSE)




