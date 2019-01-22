######################Simple bash script to run faststructure and conduct Geographic cline analysis using HZAR############
#########################################################################################################################
##author: menonm2 (menonm2@mymail.vcu.edu)
###Last modified: 18-Nov-2018
##Job submitted via VCU grid engine
##One job for all generation files per model per MC replicate
##Assuming that R is installed and is in the user's path
##########################################################################################################################


#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N ClineBi
#$ -V
#$ -j y
#$ -pe smp 10
#$ -cwd
#$ -e ./
#$ -o ./
#$ -l mem_free=40G

cd ~/path/to/CDPop/data_July2018/ModelB_ItoIII_RD1_i_Genes_1539208533/batchrun0mcrun0

######Convert CDPop genetic output for input into faststructure########
mkdir strFiles
R --no-save <  CDPoptoSTR.R


####Run faststructure##############
##NOTE: In my case faststructure is set as a seperate environment, hence the `activate faststructure` command

unset module
. activate faststructure


for f in strFiles/*.str; do
	tab=$'\t'
#remove NA in input file and replace by tab
	sed 's/NA/${tab}&/g' "$f" > strFiles/tmp
	mv strFiles/tmp "$f"
	echo "$f"
	ID="$(cut -d'.' -f1 <<< $f)"
	echo "$ID"
  #run faststructure with k as 2 for the two species
	python ~/g/src/fastStructure-master/structure.py -K 2 --input="$ID" --output="$ID" --format=str --full 
done 


######conduct geographic cline analysis for each iteration##########

R --no-save < HZAR_Bi.R

