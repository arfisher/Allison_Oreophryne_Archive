###############################
# Frog analysis script
#
# This script loads the processed, cleaned data, 
# does some analysis and saves the results
# to the results folder
###############################
## ---- setup -----
#load needed packages. make sure they are installed.
require(ggplot2) #for plotting
require(magrittr) #for piping
require(knitr) #for formatting output
require(dplyr)
require(phytools)
require(skimr)
require(rempsyc)
require(effectsize)
require(flextable)
require(corrr)
require(ggcorrplot)
require(FactoMineR)
require(factoextra)

#path to data and results 
data_path <- "../../Data/Processed_data/processeddata.rds"
results_path <- "../../Results/"

## ---- loaddata ----
dat2 <- readRDS(data_path)

## ---- functions ----
# function to paste path to output filenames

addpath <- function( filename, path=data_path ) {
    location <- paste( path, filename, sep="")
	return( location )
}


## ---- summarize ----
# create summary table of the data using skimr to use in paper
# variables, sample size, mean, standard error

sk <- skimr::skim(dat2)  # save skim object
sk <- as.data.frame(sk) # save as data.frame
head(sk)  # see the variable names

nrows <- dim(dat2)[1] # total number of rows
sk$N <- nrows - sk$n_missing  # sample size of each variable

## ---- summary.table ----
# select only the variable, N, mean, sd, and category counts

sk.table <- sk[c("skim_variable", "N", "numeric.mean", "numeric.sd")]
names(sk.table) <- c("Variable", "N", "Mean", "SE") # rename SD as SE
sk.table$SE <- sk.table$SE/sqrt(sk.table$N) # calculate SE

options(knitr.kable.NA = "")
knitr::kable(sk.table, digits=2)


# save summary table
saveRDS(sk.table, file = addpath("summary_table.rds", results_path))

## ---- summarystats ------

#lets see means for each species
dat2$gensp <- paste0(dat2$genus, dat2$species) #combine genus and species names to sort species by genus

spec.means <- dat2 %>%                                      
  group_by(gensp) %>%                         
  summarise_at(vars(SVL, finger, toe, ft, anterior, posterior, ap),          
               list(mean)) 
print(spec.means)

# and for each genus
gen.means <- dat2 %>%                                      
  group_by(genus) %>%                         
  summarise_at(vars(SVL, finger, toe, ft, anterior, posterior, ap),          
               list(mean)) 
print(gen.means)
## ---- genus_ft_boxplot ----

# plot to screen
with(dat2, boxplot(ft ~ genus))

# plot to .png file, can also do pdf using `pdf()` function 
png(filename = addpath("genus_ft_boxplot.png", results_path))
  with(dat2, boxplot(ft ~ genus))
dev.off()

## ---- genus_ap_boxplot ----

# plot to screen
with(dat2, boxplot(ap ~ genus))

# plot to .png file, can also do pdf using `pdf()` function 
png(filename = addpath("genus_ap_boxplot.png", results_path))
  with(dat2, boxplot(ap ~ genus))
dev.off()


## ---- ttests ------
# trying package rempsyc to make pretty ttest tables, grouped by genus

t.test.results <- nice_t_test(
  data = dat2,
  response = names(dat2)[4:10],
  group = "genus",
  warning = FALSE)
t.test.results

my_table <- nice_table(t.test.results)
my_table

# t-test results are really only valid for the ratios ft and ap because finger/toe size and anterior/posterior palatal groove widths likely scale to body size, should be standardized

## ---- PCA --------
#let's try a PCA, this comes from a datacamp tutorial <https://www.datacamp.com/tutorial/pca-analysis-r>


str(dat2) # show summary of data
colSums(is.na(dat2)) # check for nulls
dat3 <- dat2[,4:11] # take only numerical columns and gensp
dat3 %>% relocate(gensp, .before=SVL)
head(dat3)

dat4 <- scale(dat3) # normalize data
head(dat4)

corr_matrix <- cor(dat4) 
ggcorrplot(corr_matrix) # view a correlation matrix for the variables

data.pca <- princomp(corr_matrix) # calculate PCA analysis
summary(data.pca) # print summary of results

data.pca$loadings[, 1:2] 

fviz_eig(data.pca, addlabels = TRUE)

## ---- phylogenyexperiment ----
#I want to use one of the ways to plot either a heat map or dotTree as in this tutorial <http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html> 


dat3 <- as.data.frame(spec.means)
row.names(dat3) <- dat3$gensp #gensp must be rownames, no character columns allowed
dat4 <- dat3[,-c(1:2)]

tree_location <- "../../Data/Raw_data/asterophryinae_partitions.nex.timetree.nwk"
tree <- read.newick(tree_location)

species<-c("Oreophryneinornata",   "Oreophryneloriae" ,    "Oreophrynenotata"   , 
 "Oreophryneparkeri" ,   "Oreophrynebiroi" ,     "Oreophryneanamiatoi" 
,"Auparoparopenelopeia", "Auparoparoinsulana" ,  "Auparoparophoebe"    
, "Auparoparopicticrus" , "Auparoparomatawan"  ,  "Auparoparoezra"  ) #select included taxa

pruned.tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
write.tree(pruned.tree) #remove other taxa
plot(pruned.tree)


phylo.heatmap(pruned.tree,dat4,standardize=TRUE)

png(filename = addpath("phylo_heatmap.png", results_path))
  phylo.heatmap(pruned.tree,dat4,standardize=TRUE)
dev.off()
