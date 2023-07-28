# README.md

To reconstruct the phylogeny using IQTREE, use these four inputs: 

-  asterophryinae_20230527.phy = the sequence alignment in ../../mesquite_alignment folder
-  asterophryinae_partitions.nex = the sequence partitions file
-  asterophryinae_dates.txt  = the dates file with timings
-  "UMMZ219489_Scaphiophryne_marmorata" outgroup specified on the command line with -o flag

Flags:
-s  alignent.phy
-p  partitions.nex
-m  MFP+MERGE to run tree and find best model
    MF+MERGE  to find best model w/o running tree
-B  1000  to run tree with 1000 bootstrap replicates    
-o  "outgroup" 
--date  dates.txt

Use IQTREE2 and run the following command (IQTREE is run from the command line):

Template:

  iqtree2 -s alignment.phy -p partitions.nex -m MFP+MERGE -B 1000 --date dates.txt -o "outgroup" --date-tip 0

Filled in - run partitionfinder plus tree:

  iqtree2 -s asterophryinae_12232022.phy -p asterophryinae_partitions.nex -m MFP+MERGE -B 1000 --date asterophryinae_dates.txt -o "UMMZ219489_Scaphiophryne_marmorata" --date-tip 0

Filled in - run PartitionFinder only (not sure if we need the dates?)

  iqtree2 -s asterophryinae_12232022.phy -p asterophryinae_partitions.nex -m MF+MERGE --date asterophryinae_dates.txt -o "UMMZ219489_Scaphiophryne_marmorata" --date-tip 0

