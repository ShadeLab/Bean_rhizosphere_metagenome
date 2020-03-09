# Common bean root metagenome; 
# _Read processing and analysis_

This repository is related to the common bean metagenome processing and analysis. Samples derived from a field experiment done in 2018 at Montcalm research farm (MRF), MI. Samples were collected at 5 different growth stages from the rots of bean genotype 'Eclipse'. We used Illumina HiSeq 2x150 platform to sequence 27 samples.
 
The purpose of this work is to:
-understand the genetic potential of the bean root core taxa (Stopnisek and Shade ([bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.20.913202v2.abstract))
-understand changes in microbiome functions through plant development

## QC 
#### (using BBtols)


## Assembly

### MEGAHIT


### SPADES


##### Read normalization (BBnorm.sh)

Output:
```
 ***********   Pass 1   **********   


Settings:
threads:                12
k:                      31
deterministic:          true
toss error reads:       false
passes:                 1
bits per cell:          16
cells:                  467.58B
hashes:                 3
base min quality:       5
kmer min prob:          0.5

target depth:           160
min depth:              2
max depth:              200
min good kmers:         15
depth percentile:       64.8
ignore dupe kmers:      true
fix spikes:             false

Made hash table:        hashes = 3       mem = 870.26 GB        cells = 467.22B         used = 43.101%

Estimated unique kmers:         87820969379

Table creation time:            16153.397 seconds.
Writing interleaved.
Started output threads.
Table read time:                14153.906 seconds.      16110.82 kb/sec
Total reads in:                 1523508068      80.798% Kept
Total bases in:                 228031089493    80.810% Kept
Error reads in:                 773643916       50.780%
Error pairs in:                 485461566       63.729%
Error type 1:                   618300724       40.584%
Error type 2:                   143168532       9.397%
Error type 3:                   69759928        4.579%
Total kmers counted:            182282990216
Total unique kmer count:        91715231831
Includes forward kmers only.
The unique kmer estimate can be more accurate than the unique count, if the tables are very full.
The most accurate value is the greater of the two.

Percent unique:                 50.31%
Depth average:                  1.99    (unique kmers)
Depth median:                   1       (unique kmers)
Depth standard deviation:       5.42    (unique kmers)
Corrected depth average:        1.57    

Depth average:                  12.35   (all kmers)
Depth median:                   2       (all kmers)
Depth standard deviation:       596.17  (all kmers)
Approx. read depth median:      2.50

   ***********   Pass 2   **********   


Settings:
threads:                12
k:                      31
deterministic:          true
toss error reads:       false
passes:                 1
bits per cell:          16
cells:                  467.58B
hashes:                 3
base min quality:       5
kmer min prob:          0.5

target depth:           40
min depth:              2
max depth:              40
min good kmers:         15
depth percentile:       54.0
ignore dupe kmers:      true
fix spikes:             false

Made hash table:        hashes = 3       mem = 870.26 GB        cells = 467.22B         used = 30.506%

Estimated unique kmers:         56678862486

Table creation time:            12324.423 seconds.
Writing interleaved.
Started output threads.
Table read time:                10153.952 seconds.      18147.73 kb/sec
Total reads in:                 1230960560      93.569% Kept
Total bases in:                 184271182304    93.576% Kept
Error reads in:                 512365921       41.623%
Error pairs in:                 349851368       56.842%
Error type 1:                   361077691       29.333%
Error type 2:                   139573781       11.339%
Error type 3:                   50638201        4.114%
Total kmers counted:            147314416392
Total unique kmer count:        60849227326
Includes forward kmers only.
The unique kmer estimate can be more accurate than the unique count, if the tables are very full.
The most accurate value is the greater of the two.

Percent unique:                 41.31%
Depth average:                  2.42    (unique kmers)
Depth median:                   1       (unique kmers)
Depth standard deviation:       5.00    (unique kmers)
Corrected depth average:        2.14    

Depth average:                  11.26   (all kmers)
Depth median:                   3       (all kmers)
Depth standard deviation:       338.43  (all kmers)

Approx. read depth median:      3.75

Removing temp files.

Total time:                     52819.316 seconds.      7805.90 kb/sec
```
### IDBA-ud


## Assembly statistics 

for sample in $(<sampleID.txt)
do
stats.sh in=${sample}_megahit_assembly/final.contigs.fa format=6 out=${sample}_assembly_stats/${sample}.txt 
done



