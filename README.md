# Common bean root metagenome; 
# _Read processing and analysis_

This repository is related to the common bean metagenome processing and analysis. Samples derived from a field experiment done in 2018 at Montcalm research farm (MRF), MI. Samples were collected at 5 different growth stages from the rots of bean genotype 'Eclipse'. We used Illumina HiSeq 2x150 platform to sequence 27 samples.
 
The purpose of this work is to:
-understand the genetic potential of the bean rhizosphere core taxa (Stopnisek and Shade ([bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.20.913202v2.abstract))
-understand changes in microbiome functions through plant development


#### Data Folder
This folder contains information needed for the analysis in R. 
The MAG.* files contain information from the refined MAGs.
The MAGdepth.txt is a coverge file created by using following command from the _METAbat2_ tool
<jgi_summarize_bam_contig_depths --outputDepth MAGdepth.txt *.sort.bam>

The MAGcoverage_table.tsv was generated ussin the _concoct_ commnd:
<concoct_coverage_table.py combinedMAGS.bed *sort.bam > MAGcoverage_table.tsv>

