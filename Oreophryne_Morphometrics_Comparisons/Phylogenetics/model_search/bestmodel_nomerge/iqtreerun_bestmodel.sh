# find best model and fit tree
iqtree2 -s ../../../mesquite_alignment/oreophryne_senkenbergiana_alignment_20230715.nex  -p ../input/asterophryinae_15partitions.nex -pre senkenbergiana -m TEST -B 1000 -o "UMMZ219489_Scaphiophryne_marmorata" --date ../input/asterophryinae_dates.txt --date-tip 0 -mset JC,F81,K80,HKY,TN,TNe,TPM2,TPM2u,TPM3,TPM3u,TIM,TIMe,TIM2,TIM2e,TIM3,TIM3e,TVMe,TVM,SYM,GTR ; 
