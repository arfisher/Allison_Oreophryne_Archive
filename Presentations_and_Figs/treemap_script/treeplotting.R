# 1. Read in trees
# 2. Clean tree labels 
# 3. Trim tree down to 1 species per site  (210 sp tree plus new samples)
# 4. Plot tree with Type Oreophryne highlighted
# 5. Drop non-Oreophryne tips (or reduce to one representative per genus?)
# 6. Plot reduced trees
################################################################

require(treeio)
require(ggtree)
require(dplyr)
require(phytools)
require(ggplot2)
require(ape)
source("clean_functions.R") # contains no_()

### 1. Read in trees
beast <- read.beast("beast_asterophryinae12242022.nex")
iqtree <- read.iqtree("asterophryinae_partitions.nex.timetree.nwk")
d <- read.csv("Table1_oreotype.csv")          # tree metadata spreadsheet

d$gensp <- paste(d$genus, d$species, sep=" ")
d$tiplabid <- paste(d$gensp, d$id, sep=" ")

# type species and sisters
type <- c("Oreophryne A senckenbergiana", "Oreophryne A frontifasciata")
odds <- unlist(sapply(type, grep, d$tiplabid))
d$tipcol <- "black"
d$tipcol[odds] <- "blue"

if(dir.exists("out")!=TRUE) dir.create("out") # check if output directory out exists, if false create

### 2. Clean up the tip labels - extra "_" in some names

	# iqtree
namesiqtree <- tips <- get_taxa_name(ggtree(iqtree)) 
#tips <- sub("FK_6919", "BPBM16806", tips)
tips <- sub("_.*$", "", tips)

iqtips <- data.frame("id"=tips, "treelabel"=namesiqtree)
iqd <- merge(iqtips, d, by="id", all=T)  # check that dataframe matches tree

iqtree <- rename_taxa(iqtree, iqd, treelabel, label) # change tip labels from the original labels to cleaned labels

	# beast
namesbeast <- tips <- get_taxa_name(ggtree(beast))
tips <- sub("FK_6919", "BPBM16806", tips)
tips  <- sub("_.*$", "", tips)

beasttips <- data.frame("id"=tips, "treelabel"=namesbeast)
beastd <- merge(beasttips, d, by="id", all=T)  # check that dataframe matches tree

tree <- beast <- rename_taxa(beast, beastd, treelabel, label) # change tip labels from the original labels to cleaned labels
tree_names <- get_taxa_name(ggtree(tree))

### 3. Trim tree down to 1 sample per species (210 sp tree)

tree <- full_join(tree, d, by="label") # full tree with metadata
tibb <- as_tibble(tree)
ss <- distinct(tibb, genus, species, .keep_all=T)  # distinct species per site
todrop <- as_tibble(tree)$label %w/o% ss$label # drop 30 taxa at multiple sites, outgroups
todrop <- grep("fronti", todrop, invert=T, value=T )  # leave in all the O. frontifasciata, drop 28 taxa

outgroups <- c("UMMZ211174", "UMMZ211181", "UMMZ219489")
todrop_og <- grep(paste(outgroups, collapse="|"), tree_names, value=T) 	# find index to drop - low sequence info

treefull <- tree   # 240 tips
tree <- treeio::drop.tip(tree, c(todrop, todrop_og))  # ingroup tree with only one sample per species - 212 tips
write.beast(tree, file="out/beast_212_tree_asterophryinae_oreophryne.nex") #ingroup tree, 1 species per site


### 4. Plot tree with Type Oreophryne highlighted

tdat <- tree %>% as_tibble %>% as.data.frame
#tdat %>% filter(!is.na(tiplabid)) #remove rows with NA tiplabid

genera <- unique(d$genus) %w/o% c("Oxydactyla", "Sphenophryne", "Aphantophryne A", "Aphantophryne B", "Dyscophis", "Platypelis", "Scaphiophryne")   # 17 recognized genera

gc <- read.csv("gencolorABC.csv")    ## color code highlight boxes by genus
gcol <- gc$col
names(gcol) <- gc$gen         ## gcol holds the colors=values for the genera=names  

# plot no ID tree: (use label=tiplab2id for labels with IDs)
#		gensp = genus species
# 		tiplabid = genus species id
#		tiplab2	= genus species-site if at multiple sites
# 		tiplab2id = tiplab2 plus id

p <- ggtree(tree) +     # 212 taxa
	 geom_tiplab(aes(label=gensp, color=tipcol), 
	 					size=1.2, 
	 					offset=.05, 
	 					fontface='italic' ) + 
	 scale_color_manual(values = c("black", "blue")) +
     scale_x_continuous(limits = c(0, 20), 
     					breaks=c(-0.71, 4.71, 9.71, 14.71, 19.71),  
     					labels = c(20,15,10,5,0)) +
     coord_cartesian(clip="off") + 
	 theme_tree2(legend.position="none") 

pdf(file="out/FullTree.pdf", height=10, width=3)
  print(p)
dev.off()


### Identifying genera

# find basal node for each genus to use in highlighting clades
# first make tipvec, which excludes the oddball species from defining genera (probably misindentifications or non-monophyletic)
oddballs <- c("Paedophryne dekot", "Aphantophryne", "Copiula tyleri", "Oreophryne A graminis")
tipvec <- grep(paste(c(nospaces(oddballs),"Copiula_sp.7"), collapse="|"), get_taxa_name(ggtree(tree)), value=T, invert=T) 

# gennode = vector containing node numbers at root of each genus (excluding oddballs)
gennode <- sapply(nospaces(genera), function(x)   # for each set of tips per genus, find MRCA
                            MRCA(tree,   
                              grep(x, 
								tipvec,
								value=T)
							)
		         )

# test tree if you want to show node numbers
#test <- p + geom_text2(aes(label=node), hjust=-.2, size=1.5) 
     
## adds genus highlight boxes  
q <- p + geom_hilight(node=gennode, fill=gcol[no_(names(gennode))], alpha=0.25)

pdf(file="out/FullTreeGenera.pdf", height=10, width=4)
  print(q)
dev.off()

### highlight just Oreophryne A, B 

justoreo <- grep("Oreophryne", names(gennode), value=T )

qo <- p + geom_hilight(node=gennode[justoreo], fill=gcol[no_(names(gennode[justoreo]))], alpha=0.25) 

pdf(file="out/FullTreeOreoHighlighted.pdf", height=10, width=4)
  print(qo)
dev.off()

# with support values

pdf(file="out/FullTreeOreoHighlightedSupp.pdf", height=10, width=4)
  print(qo + geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=2))
dev.off()


### 5. Plot trees, dropping non-Oreophryne tips (or reduce to one representative per genus?)
### This code doesnʻt modify the tree object, just drops it from the plot

# Collapse Paedophryne 

p <- ggtree(tree) +     # 212 taxa
	 geom_tiplab(aes(label=gensp, color=tipcol), 
	 					size=3, 
	 					offset=.05, 
	 					fontface='italic' ) + 
	 scale_color_manual(values = c("black", "blue")) +
     scale_x_continuous(limits = c(0, 50), 
     					breaks=c(-0.71, 4.71, 9.71, 14.71, 19.71),  
     					labels = c(20,15,10,5,0)) +
     coord_cartesian(clip="off") + 
	 theme_tree2(legend.position="none") 

# This stores the node that is the MRCA of all the rest of the tree that we are dropping. If the topology changes, what is collapsed will need to change
restoftree <- MRCA(tree, gennode[c("Choerophryne", "Hylophorbus")])

r <- p %>% collapse(gennode["Paedophryne"]) %>% collapse(gennode["Cophixalus"]) %>% collapse(restoftree) +
  geom_point2(aes(subset=(node==gennode[c("Paedophryne")])), shape=21, size=5, fill=gcol[c("Paedophryne")], alpha=0.5) +
  geom_point2(aes(subset=(node==gennode[c("Cophixalus")])), shape=21, size=5, fill=gcol[c("Cophixalus")], alpha=0.5) +
  geom_point2(aes(subset=(node==restoftree)), shape=21, size=5, fill=gcol[c("anc")], alpha=0.5) + 
  geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3)  # add the support values

pdf(file="out/OreoTreeSupp.pdf", height=10, width=4)
  print(r)
dev.off()


### 6. Drop non-Oreophryne tips from tree, plot map 

tree <- iqtree		# use iqtree for now, can change to beast tree later
#dat <- d[grep("Oreophryne A|Oreophryne B", d$genus),]   # metadata containing only Oreoprhyne species. Use to drop tips and get gps
dat <- read.csv("Table1_onlyoreo.csv")   

tips <- get_taxa_name(ggtree(tree))
ids  <- sub("[_-].*$", "", tips)
tipdat <- data.frame (labels=tips, id=ids)

tree <- rename_taxa(tree, tipdat, labels, id)  # tip labels now only have id

tokeep <- dat$id
todrop <- ids %w/o% dat$id
tree <- treeio::drop.tip(tree, todrop)
treep <- as.phylo(tree)   # phytools requires trees in phylo format

rownames(dat) <- dat$id # paste id numbers into the rownames of the gps data

obj<-phylo.to.map(treep,dat[c("latitude", "longitude")],plot=FALSE) # the object that contains the map+phylogeny
plot(obj,type="phylogram",asp=0.75,mar=c(0.1,0.5,3.1,0.1), xlim=c(125,155), ylim=c(-15,3), fsize=.5)

#makepdf
pdf(file="out/Oreomap.pdf", height=6, width=10)
  plot(obj,type="phylogram",asp=1,mar=c(0.1,0.5,2.1,0.1), tree.mar=c(0.1,0.1,0.5,1.1), xlim=c(120,155), ylim=c(-15,3), fsize=.5, split = c(.25,.75))
dev.off()

### 7. Plot maps + phylogenies for the two clades Oreophryne A (Oreophryne) and Oreophryne B (Auparoparo)

## just Oreophryne 
oreoA <- dat[dat$genus == "Oreophryne A", ]   # subset metadata to only oreoA
treeA <- keep.tip(treep, oreoA$id)                     # drop tips keeping only oreoA

objA <- phylo.to.map(treeA, oreoA[c("latitude", "longitude")], plot=FALSE)
plot(objA,type="phylogram",asp=1,mar=c(0.1,0.5,2.1,0.1), tree.mar=c(0.1,0.1,0.5,1.1), xlim=c(120,155), ylim=c(-15,3), fsize=.5, colors=gcol["Oreophryne A"], split = c(.25,.75))

pdf(file="out/OreoAmap.pdf", height=6, width=10)
  plot(objA,type="phylogram",asp=1,mar=c(0.1,0.5,2.1,0.1), tree.mar=c(0.1,0.1,0.5,1.1), xlim=c(120,155), ylim=c(-15,3), fsize=.5, colors=gcol["Oreophryne A"], split = c(.25,.75))
dev.off()
#for green map use colors of "Hylophorbus"

## just Auparoparo 
oreoB <- dat[dat$genus == "Oreophryne B", ]
treeB <- keep.tip(treep, oreoB$id)

objB <- phylo.to.map(treeB, oreoB[c("latitude", "longitude")], plot=FALSE)
plot(objB,type="phylogram",asp=1,mar=c(0.1,0.5,2.1,0.1), tree.mar=c(0.1,0.1,0.5,1.1), xlim=c(120,155), ylim=c(-15,3), fsize=.5, split = c(.25,.75))


pdf(file="out/OreoBmap.pdf", height=6, width=10)
  plot(objB,type="phylogram",asp=1,mar=c(0.1,0.5,2.1,0.1), tree.mar=c(0.1,0.1,0.5,1.1), xlim=c(120,155), ylim=c(-15,3), fsize=.5, colors=gcol["Oreophryne B"], split = c(.25,.75))
dev.off()

#for red map remove colors command (default is red)