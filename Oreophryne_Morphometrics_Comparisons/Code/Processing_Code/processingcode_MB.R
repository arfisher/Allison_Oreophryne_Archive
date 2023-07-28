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

## ---- functions -----

# function to paste path to output filenames
addpath <- function( filename, path=data_path ) {
    location <- paste( path, filename, sep="")
	return( location )
}

# function to create a new variable from multiple dummy variable columns
collapsevars <- function( dat) {
	newvar <- names(dat)[max.col(dat)]
	return(newvar)
}

# x without y function 
 "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y

## ---- loaddata --------
data_location <- "../../Data/Raw_Data/Oreophryne_Character_Measurements_Raw.csv"
data_path <- "../../Data/Raw_Data/"
results_path <- "../../Results/"

rawdata <- read.csv(data_location, check.names=FALSE)
head(rawdata)

# adding in score data
scores <- read.csv(paste(data_path, "Oreophryne_character_scores.csv", sep=""))
print(scores)
names(scores) <- gsub(" ", "_", names(scores))
names(scores) <- gsub("_call", "", names(scores))
names(scores)[names(scores)=="SVL"] <- "svl"
names(scores)[names(scores)=="X3Fpad.4Tpad"] <- "fpadtpad"

# collapse multiple column variables to single column
pecvars <- c("cartilaginous", "ligamentous")
callvars <- c("honk", "whinny", "rattle", "peeping", "other")

scores$pectoral_type <- collapsevars( scores[pecvars] )
scores$call_type <- collapsevars( scores[callvars] )

#scores <- scores[names(scores) %w/o% c(pecvars, callvars)]

scores <- scores %>% 
				mutate(toe_length = if_else( toe_length==3, TRUE, FALSE )) %>%
				mutate(pupil_horizontal = if_else( pupil_horizontal==3, TRUE, FALSE )) %>%
				mutate(teeth = if_else( teeth==3, TRUE, FALSE )) %>%
				mutate(cartilaginous = if_else( cartilaginous ==3, TRUE, FALSE )) %>%
				mutate(ligamentous = if_else( ligamentous ==3, TRUE, FALSE )) %>%
				mutate(honk = if_else( honk ==3, TRUE, FALSE )) %>%
				mutate(whinny = if_else( whinny ==3, TRUE, FALSE )) %>%
				mutate(rattle = if_else( rattle ==3, TRUE, FALSE )) %>%
				mutate(peeping = if_else( peeping ==3, TRUE, FALSE )) %>%
				mutate(other = if_else( other ==3, TRUE, FALSE ))

# summarize by genus and species
mscore <- scores %>% group_by(genus, species) %>% 
		summarise(across(svl:other, ~ mean(.x, na.rm = TRUE)))

# summarize by genus
gscore <- mscore %>% group_by(genus) %>% 
		summarise(across(svl:other, ~ mean(.x, na.rm = TRUE)))

# plot calls, organize data by call type

# likert plot
calls <- gscore[c("genus", callvars)]
likert(species ~.|genus, layout=c(1,2), calls, positive.order = TRUE, 
       scales=list(y=list(relation="free")),
       strip.left=strip.custom(bg="gray97"),
       strip=FALSE,
       as.percent = "noRightAxis", ReferenceZero = 2.5,
       main = 'Call types', 
       ylab = "Genus", xlab = "Percentage",
       sub= list("Angry Level Rating",x=unit(.6, "npc")))

# stacked bar plot
gcall <- mscore %>% group_by(genus) %>% 
		summarise(across(honk:other, ~ sum(.x, na.rm = TRUE))) %>%
		pivot_longer(-genus, names_to="call_type", values_to="n")

ggplot(gcall, aes(fill=call_type, y=n, x=genus)) + 
    geom_bar(position="fill", stat="identity")


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
colnames(dat) <- c("genus","species", "BPBM", "svl", "finger", "toe", "ft", "anterior", "posterior", "ap")

head(dat)
## ---- merge ----
#merge measurement with score data
dat2 <- merge(dat, scores, by = c('BPBM', 'genus', 'species', 'svl'), all=TRUE)
dat2$gensp <- with(dat2, paste(genus, species, sep="_"))

## ---- exploredata --------

# look at the data
skimr::skim(dat2)

## ---- exploratoryplots --------
# create scatter plots of SVL vs finger:toe ratio and SVL vs anterior width:posterior width of palatal groove, colored by genus
svl.vs.ft <- ggplot(data = dat2) + geom_point(aes(x = svl, y = ft, col=dat2$genus))
svl.vs.ft

png(filename = addpath("svl_ft.png", results_path))
  svl.vs.ft
dev.off()

svl.vs.ap <- ggplot(data = dat2) + geom_point(aes(x = svl, y = ap, col=dat2$genus))
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