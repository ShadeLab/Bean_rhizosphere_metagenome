# Common bean root metagenome; 
#_Read processing and analysis_

This repository is related to the common bean metagenome processing and analysis. Samples derived from a field experiment done in 2018 at Montcalm research farm (MRF), MI. Samples were collected at 5 different growth stages from the rots of bean genotype 'Eclipse'. We used Illumina HiSeq 2x150 platform to sequence 27 samples.
 
The purpose of this work is to:
-understand the genetic potential of the bean root core taxa (Stopnisek and Shade ([bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.20.913202v2.abstract))
-understand changes in microbiome functions through plant development

##QC 
####(using BBtols)


## Assembly

### MEGAHIT


### SPADES


## Assembly statistics 

for sample in $(<sampleID.txt)
do
stats.sh in=${sample}_megahit_assembly/final.contigs.fa format=6 out=${sample}_assembly_stats/${sample}.txt 
done

