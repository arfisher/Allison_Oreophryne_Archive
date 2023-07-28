###############################
# Frog processing script
#
#this script loads the raw data, processes and cleans it 
#and saves it as Rds file in the Processed_data folder

## ---- packages --------
#load needed packages. make sure they are installed.
require(ggplot2) #exploring plots
require(dplyr) #for data processing/cleaning
require(skimr) #for nice visualization of data 

# function to paste path to output filenames

addpath <- function( filename, path=data_path ) {
    location <- paste( path, filename, sep="")
	return( location )
}


## ---- loaddata --------
data_location <- "../../Data/Raw_data/Oreophryne_Character_Measurements.csv"
data_path <- "../../Data/Raw_data/"
results_path <- "../../Results/"

rawdata <- read.csv(data_location, check.names=FALSE)
head(rawdata)

# adding in score data
scores <- read.csv(paste(data_path, "Oreophryne_character_scores.csv", sep=""))
print(scores)
names(scores) <- gsub(" ", "_", names(scores))

# view data dictionary
dictionary <- read.csv(paste(data_path, "datadictionary.csv", sep=""))
print(dictionary)

## ---- cleanup --------
head(rawdata) #notice because of the format of ImageJ output and the way I converted pixels to mm, there are 5 rows with NA for each specimen

#Making a new rectangular dataframe for only converted measurements (mm not pixels) for each specimen

# exclude rows with NAs
nas <- which( is.na(rawdata$BPBM) ) # find which rows have NA 
dat <- rawdata[-nas,] # exclude these rows
head(dat)

#removing unnecessary columns, artifacts from importing pixel measurements from ImageJ
dat <- dat[-c(5:13)]

#cleaning up column names
colnames(dat) <- c("genus","species", "BPBM", "SVL", "finger", "toe", "ft", "anterior", "posterior", "ap")

head(dat)
## ---- merge ----
#merge measurement with score data
dat2 <- merge(dat, scores, by = c('BPBM', 'genus', 'species', 'SVL'), all=TRUE)

## ---- exploredata --------

# look at the data
skimr::skim(dat2)

## ---- exploratoryplots --------
# create scatter plots of SVL vs finger:toe ratio and SVL vs anterior width:posterior width of palatal groove, colored by genus
svl.vs.ft <- ggplot(data = dat2) + geom_point(aes(x = SVL, y = ft, col=dat2$genus))
svl.vs.ft

png(filename = addpath("svl_ft.png", results_path))
  svl.vs.ft
dev.off()

svl.vs.ap <- ggplot(data = dat2) + geom_point(aes(x = SVL, y = ap, col=dat2$genus))
svl.vs.ap

png(filename = addpath("svl_ap.png", results_path))
  svl.vs.ap
dev.off()

# just for fun... ft vs ap by genus
ft.vs.ap <- ggplot(data = dat2) + geom_point(aes(x = ft, y = ap, col=dat2$genus))
ft.vs.ap

png(filename = addpath("ft_ap.png", results_path))
  ft.vs.ap
dev.off()

#density plots
ft.dens <- dat2 %>%    # CC by species
        ggplot( aes(x=`ft`)) + 
		geom_density( aes(fill=genus), alpha=.5)
ft.dens

png(filename = addpath("ft_dens.png", results_path))
  ft.dens
dev.off()

ap.dens <- dat2 %>%    # CC by species
        ggplot( aes(x=`ap`)) + 
		geom_density( aes(fill=genus), alpha=.5)
ap.dens

png(filename = addpath("ap_dens.png", results_path))
  ap.dens
dev.off()
# we can see some trends but no clear distinctions yet

## ---- savedata --------
processeddata <- dat2      # change if you did more steps

# location to save file
save_data_location <- "../../Data/Processed_data/processeddata2.rds"
saveRDS(processeddata, file = save_data_location)

save_data_location_csv <- "../../Data/Processed_data/processeddata2.csv"
write.csv(processeddata, file = save_data_location_csv, row.names=FALSE)

print(dat2 )