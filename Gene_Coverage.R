library(tidyverse)
library(vegan)
library(gplots)

setwd('Documents/git/Bean_rhizosphere_metagenome/')

geneList=read.table('../LargeFiles/BeanMetaG/functionsAnvio.txt', header=T, sep='\t')
head(geneList)
geneCov=read.table('../LargeFiles/BeanMetaG/AnvioCoverage-GENE-COVERAGES.txt', sep='\t', header=T)
head(geneCov)

mapFile=data.frame(sampleID=colnames(geneCov)[-1], time=NA)
mapFile=mapFile %>%
  separate(sampleID,sep='_', into = c('a', 'b', 'growthStage', 'plot', 'plant',
                                      'c', 'd', 'e'), remove = F)
mapFile=mapFile[-c(2,3,7,8,9,10)]
head(mapFile)

#' Subset the geneList table to the levels of description
#' [-2] removes the source column
functionList=geneList[geneList$source == "COG20_FUNCTION",][-2]
categoryList=geneList[geneList$source == "COG20_CATEGORY",][-2]
pathwayList=geneList[geneList$source == "COG20_PATHWAY",][-2]

#' COG CATEGORY
#' Subset the gene coverage table to genes predicted and annotated with category
CATgene=geneCov[geneCov$key %in% categoryList$gene_callers_id,]

#' Aggregate the genes by category
#' Using melt() function to make a long table first and aggregate
CATgeneLong=CATgene %>%
  left_join(categoryList[c(1,2)], by=c( 'key' = 'gene_callers_id')) %>%
  pivot_longer(
    cols=QTRIMMED_MRF_FLOW_1_1_S17_L007_SORT:QTRIMMED_MRF_V5_3_2_S34_L008_SORT,
    values_to='coverage',
    names_to='sampleID')

head(CATgeneLong)

categoryDF=CATgeneLong %>%
  group_by(accession, sampleID) %>%
  summarise(coverageSum=sum(coverage))

#' How many unique categories?
length(unique(categoryDF$accession))
# 1744

#' *****************************************************************************
#' COG FUNCTION
#' Subset the gene coverage table to genes predicted and annotated with function

FUNgene=geneCov[geneCov$key %in% functionList$gene_callers_id,]

#' Aggregate the genes by category
#' Using melt() function to make a long table first and aggregate
FUNgeneLong=FUNgene %>%
  left_join(functionList[c(1,2)], by=c( 'key' = 'gene_callers_id')) %>%
  pivot_longer(
    cols=QTRIMMED_MRF_FLOW_1_1_S17_L007_SORT:QTRIMMED_MRF_V5_3_2_S34_L008_SORT,
    values_to='coverage',
    names_to='sampleID')

head(FUNgeneLong)

functionDF=FUNgeneLong %>%
  group_by(accession, sampleID) %>%
  summarise(coverageSum=sum(coverage))

#' How many unique functions?
length(unique(functionDF$accession))
#' 10772
#' 
#' #' *****************************************************************************
#' COG PATHWAYS
#' Subset the gene coverage table to genes predicted and annotated with function

PATHgene=geneCov[geneCov$key %in% pathwayList$gene_callers_id,]

#' Aggregate the genes by category
#' Using melt() function to make a long table first and aggregate
PATHgeneLong=PATHgene %>%
  left_join(pathwayList[c(1,2)], by=c( 'key' = 'gene_callers_id')) %>%
  pivot_longer(
    cols=QTRIMMED_MRF_FLOW_1_1_S17_L007_SORT:QTRIMMED_MRF_V5_3_2_S34_L008_SORT,
    values_to='coverage',
    names_to='sampleID') %>%
  left_join(mapFile[c(1,2,3,4)]) %>%
  mutate(sampleID=paste(growthStage, plot, plant, sep='.'))

PATHgeneLong$accession=vapply(strsplit(PATHgeneLong$accession,"!!!"), `[`, 1, FUN.VALUE=character(1))
head(PATHgeneLong)

pathwayDF=PATHgeneLong %>%
  group_by(accession, sampleID) %>%
  summarise(coverageSum=sum(coverage))

#' How many unique functions?
length(unique(pathwayDF$accession))
#' 947

PATHwide=pathwayDF %>%
  pivot_wider(names_from=sampleID, values_from=coverageSum)

PATH.matrix=data.matrix(PATHwide, rownames.force = NULL)
rownames(PATH.matrix) = PATHwide$accession
PATHmatrix=PATH.matrix[,-1]
head(PATH.matrix)

norm.PATH=decostand(PATHmatrix, method='standardize', MARGIN = 1)

hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

heatmap.2(norm.PATH, col=hc(100), trace='none', ,density.info="none", 
          sepcolor="black", Colv = F)
