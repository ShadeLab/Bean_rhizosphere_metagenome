# Working with trees
library(ape)
library(phytools)
# Plotting Grammar of Graphics
library(ggplot2)
# Tree visualization
library(ggtree)
# String maniputation
library(stringr)
# Data structure manipulation
library(reshape2)
# Multi-panel figure layouts
library(cowplot)
library(grid)
# Color maniputation
library(RColorBrewer)
# Pathway analysis with KEGG Orthologs
library(pathview)
#plot circular graphs and contig maps
library(circlize)
# Suite of packages for data manipulation and visualization
library(tidyverse)
library(treeio)
library(gameofthrones)
library(ggtreeExtra)

################################################################################
# Downloading custom functions
################################################################################

setwd("~/Documents/git/Bean_rhizosphere_metagenome/")
tree.path <- "data/Metabat1500.gtdbtk.bac120.classify.tree"

#Use read.tree(<tree.path>) to load tree file from path.
#Important: Make sure you use the "ape" version of the read.tree() function
tree <- ape::read.tree(tree.path)

#Inspect tree object
tree

core16s <- read.delim("data/core16s.txt", na.strings = '')
core16s <- core16s %>%
  mutate(Kingdom=paste("k__",Kingdom,sep=''),
         Phylum=paste("p__",Phylum,sep=''),
         Class=paste("c__",Class,sep=''),
         Order=if_else(is.na(Order), Order, false=paste("o__",Order,sep='')),
         Family=if_else(is.na(Family), Family, false=paste("f__",Family,sep='')),
         Genus=if_else(is.na(Genus), Genus, false=paste("g__",Genus,sep='')),
         Species=if_else(is.na(Species), Species, false=paste("s__",Species,sep='')))

################################################################################
# Plot whole tree with ggtree
################################################################################

GTDB.taxonomy.path <- "data/gtdb_taxonomy.tsv"
metabat.taxonomy.path <- "data/metabat1500.gtdbtk.classification.csv"
#Load data
GTDB.taxonomy <- read.csv(GTDB.taxonomy.path, sep="\t", stringsAsFactors = FALSE, head=F)
names(GTDB.taxonomy)[1]='GenomeID'
names(GTDB.taxonomy)[2]='Taxonomy'
GTDB.taxonomy=separate(data=GTDB.taxonomy, col = Taxonomy, into = c('Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species'), sep = ";")

metabat.taxonomy <- read.csv(metabat.taxonomy.path, stringsAsFactors = FALSE)
names(metabat.taxonomy)[1]='GenomeID'
names(metabat.taxonomy)[2]='Taxonomy'
metabat.taxonomy=separate(data=metabat.taxonomy, col = Taxonomy, into = c('Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species'), sep = ";")

#Inspect data
names(GTDB.taxonomy)
names(metabat.taxonomy)

#Merge datasets
taxonomy.all <- rbind(GTDB.taxonomy,metabat.taxonomy)
taxonomy.all<- taxonomy.all[taxonomy.all$Kingdom !='d__Archaea',] #removing Archaea from the taxonomy file

#Isolate and inspect tip labels. Each tip is a genomic placement. 
tips <- tree$tip.label
head(tips)

#Filter only taxa present in the tree
filtered.taxonomy=taxonomy.all[taxonomy.all$GenomeID %in% tips,]
tips.keep <- tips[tips %in% filtered.taxonomy$GenomeID]
tree.bact <- ape::keep.tip(tree, tips.keep)
tree <- tree.bact

#One more thing: adding bootstrap values to trees
#tree$node.label <- parse_bootstraps(tree, method = "parse")

#Find and parse tree boostraps
#bs_count <- parse_bootstraps(tree, method = "count")
#tail(bs_count)

Class.nodes <- collapse_nodes("Class", tree, filtered.taxonomy)
coreClass_node <- Class.nodes[which(Class.nodes$Group %in% unique(core16s$Class)),]$Node

Phyla.nodes <- collapse_nodes("Phylum", tree, filtered.taxonomy)

acido_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Acidobacteriota"),]$Node
actiono_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Actinobacteriota"),]$Node
bact_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Bacteroidota"),]$Node
chloro_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Chloroflexota"),]$Node
plan_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Planctomycetota"),]$Node
proteo_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Proteobacteria"),]$Node
ver_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Verrucomicrobiota"),]$Node

sp.nodes <- collapse_nodes("Species", tree, filtered.taxonomy)
markingData <- tree$tip.label
markingData=data.frame(genomes=markingData, group=NA) %>%
  mutate(groups=ifelse(genomes %in% metabat.taxonomy$GenomeID, 'MAG','other'))

row.names(markingData) <- NULL
#Plot tree with bootstraps as nodepoint

pal <- got(7, option = "Daenerys") # generating 7 colors for highlighting phyla

tree.plot <- ggtree(tree, layout = 'circular', color="darkgrey", size=.2) + 
  xlim(0,4)  +
  geom_cladelabel(node=acido_node, label="Acidobacteriota",align=TRUE, offset=.5) +
  geom_cladelabel(node=actiono_node, label="Actinobacteriota",align=TRUE, offset=.5) +
  geom_cladelabel(node=bact_node, label="Bacteroidota",align=TRUE, offset=.5) +
  geom_cladelabel(node=chloro_node, label="Chloroflexota",align=TRUE, offset=.5) +
  geom_cladelabel(node=plan_node, label="Planctomycetota",align=TRUE, offset=.5) +
  geom_cladelabel(node=proteo_node, label="Proteobacteria",align=TRUE, offset=.5) +
  geom_cladelabel(node=ver_node, label="Verrucomicrobiota",align=TRUE, offset=.5) +
  geom_hilight(node=acido_node, fill=pal[1], alpha=.5) +
  geom_hilight(node=actiono_node, fill=pal[2], alpha=.5) +
  geom_hilight(node=bact_node, fill=pal[3], alpha=.5) +
  geom_hilight(node=chloro_node, fill=pal[4], alpha=.5) +
  geom_hilight(node=plan_node, fill=pal[5], alpha=.5) +
  geom_hilight(node=proteo_node, fill=pal[6], alpha=.5) +
  geom_hilight(node=ver_node, fill=pal[7], alpha=.5) 

treeFig <- tree.plot %<+% markingData + 
  geom_tippoint(aes(fill=groups, color=groups), size=1) +
  scale_fill_manual(values = c('gold', 'transparent')) +
  scale_color_manual(values = c('gold', 'transparent')) 

#' First figure showing all the MAGs binned by metabat and 1500bp threshold
treeFig


# Now plot only the MAGs without the references and include the completeness bars
# Import checkM results
metabat.checkM <- read.csv('data/metabat1500_checkm.csv', header = T)
head(metabat.checkM)
class(metabat.taxonomy$Kingdom)
dat1 <- left_join(metabat.checkM, metabat.taxonomy, by="GenomeID")
dat1 <- dat1[complete.cases(dat1$Kingdom), ]
dat1$GenomeID=as.character(dat1$GenomeID)
# Subset tree to MAGs only
tips.mags <- tips[tips %in% dat1$GenomeID]
tree.metabat <- ape::keep.tip(tree, tips.mags)

dat2  <- melt(dat1[c(1:4,6)], id.vars = c("GenomeID",'Phylum'))
head(dat2)

# Make same order of tip labels and GenomeIDs
dat2$GenomeID <- factor(dat2$GenomeID, 
                                   levels = tree.metabat$tip.label)
dat2$variable <- factor(dat2$variable, 
                                   levels = rev(unique(dat2$variable)))

completness <- dat2[dat2$variable == 'Completeness',]
contamination <- dat2[dat2$variable == 'Contamination',]
names(completness)[4]='completness'
names(contamination)[4]='contamination'

df <- left_join(completness[-3], contamination[-3])
head(df)

completnessFig <- ggplot(df, aes(x = fct_reorder(GenomeID, -completness), 
                                 y = completness), fill='cadetblue') +
  geom_bar(stat='identity', color = "black") + #The "black" color provides the border. 
  theme_classic()+
  geom_hline(yintercept=70, color='grey70', linetype='dashed')+
    theme(legend.position = "none",
          axis.text.y= element_blank()) + #We don't need a legend for these data
  coord_flip() +
  ylab("Completeness")+
  xlab('GenomeID') #Add label 

library(scales)

contaminationFig <- ggplot(df, aes(x =fct_reorder(GenomeID, -completness), 
                                   y = contamination)) +
  geom_bar(stat = "identity", color = "black", fill='brown3') + #The "black" color provides the border. 
  theme_classic()+
  theme(legend.position = "none", 
        axis.text.y= element_blank(),
        axis.title.y = element_blank()) + #We don't need a legend for these data
  geom_hline(yintercept=10, color='grey70', linetype='dashed')+
  coord_flip() +
  scale_y_continuous(limits=c(0,100),oob = rescale_none)+
  ylab("Contamination")+
  xlab('GenomeID') #Add label 

grid.arrange(completnessFig, contaminationFig, widths=c(2,2))

#Where are the thresholds
sum(df$completness>70)  # 66 MAGs with 70% completeness or higher
sum(df$contamination<10) # 89 MAGs with 10% contamination or lower

sum(df$completness>70 & df$contamination<10) # 18 MAGs with 70%> completeness and 10%< contamination
sum(df$completness>50 & df$contamination<10) # 29 MAGs with 70%> completeness and 10%< contamination

Phyla.nodes <- collapse_nodes("Phylum", tree.metabat, dat1)
unique(dat1$Phylum)
metabatTree <- ggtree(tree.metabat, ladderize=F) + geom_tiplab(size = 3)
Fig2 <- clade_labels(metabatTree, Phyla.nodes, tiplimit = 1, fontsize = 4)
  

ggtree(tree.metabat, layout= 'circular', branch.length='none') +
  geom_fruit(
    data=completness, # The abundance of dat1 will be mapped to x, 
    geom=geom_bar,
    mapping=aes(y=GenomeID, x=Abundance, fill=Phylum),
    stat="identity")
  
  geom_fruit(data=completness, 
             geom=geom_bar,
             mapping=aes(fill=value),
    starstroke=0.2
  )
  
  
  

################################################################################
# Plot the reduced tree with ggtree
################################################################################

#Plot the tree using rectangular layout
tree.plot <- ggtree(tree + xlim(0,4)
tree.plot

#Plot the tree using ciculat layout
ggtree(tree.subset, layout="circular") 
  

#add tip labels with geom_tiplab()
tree.plot + geom_tiplab(size = 2)

################################################################################
# Add phyla labels with geom_cladelabel()
################################################################################

#Assign our phylogenetic clade of interest to a variable called "Clade"
Clade <- unique(core16s$Phylum)

#Here we combine a few base R functions to subset the GTDB taxonomy table to only: 
GTDB.taxonomy.clade <- filteredGTDB.taxonomy[which(filteredGTDB.taxonomy[,"Phylum"] %in% Clade),]
corePhyla_tips <- GTDB.taxonomy.clade$GenomeID

#find the common ancestor nodes for a list of tree tips. 
corePhyla_node <- findMRCA(whole_tree, corePhyla_tips, type="node")


#Use custom function: collapse_nodes() to apply findMRCA() over the whole tree to collect nodes for all phyla.
Phyla.nodes <- collapse_nodes("Phylum", whole_tree, filteredGTDB.taxonomy)

#Almost done, next we can highlight clades of interest. 
corePhylum_node <- Phyla.nodes[which(Phyla.nodes$Group %in% Clade), ]$Node

################################################################################
# Plot the tree with ggtree
################################################################################

#We can run geom_balance() to produce colored ranges for all phyla in our dataframe in a loop. 

colorS <- randomColor(length(unique(Phyla.nodes$Group)))

Phyla.nodes <- collapse_nodes("Phylum", whole_tree, filteredGTDB.taxonomy)

for(i in 1:nrow(Phyla.nodes)) {
  if(Phyla.nodes$TipCount[i] > 10) {
    tree.plot <- tree.plot + 
      geom_balance(node=Phyla.nodes$Node[i], fill=colorS[i], color='white', alpha=0.6, extend=0.6)  +
      geom_cladelabel(node = Phyla.nodes$Node[i], Phyla.nodes$Group[i], offset = 0.6, fontsize = 2)
  }
}
tree.plot


#Can we do this for a circular plot?
tree.plot <- ggtree(whole_tree, layout="circular")
tree.plot <- clade_colors(whole_tree, Phyla.nodes, colorS)
tree.plot

#(Looks pretty badass to me!)

#One more thing: adding bootstrap values to trees
whole_tree$node.label <- parse_bootstraps(whole_tree, method = "parse")

#Find and parse tree boostraps
bs_count <- parse_bootstraps(whole_tree, method = "count")

#Plot tree with bootstraps as text
ggtree(whole_tree,layout="circular") + xlim(0,3) + geom_text(aes(label = bs_count), size = 2)

#Plot tree with bootstraps as nodepoint
tree.plot <- ggtree(whole_tree) + 
  xlim(0,4) + 
  geom_nodepoint(aes(subset= bs_count >= 90), 
                 fill = "cadetblue",
                 size=1.5, 
                 alpha = 0.5, 
                 shape = 21)

tree.plot


# And add phyla labels:
tree.plot <- clade_labels(tree.plot, Phyla.nodes, tiplimit = 5, fontsize = 2)
tree.plot


#add a scale representing the substitution distance

tree.plot + geom_treescale(x=0, y=length(tree.subset$tip.label)-50, width=0.2, offset = 10)

################################################################################
# Plot the reduced tree with only refined MAGs with ggtree
################################################################################
#Load data
mag.taxonomy <- read.csv("data/MAG.gtdbtk.classification.csv", stringsAsFactors = FALSE, na.strings = '')
mag.taxonomy <- mag.taxonomy %>%
  mutate(Kingdom=paste("k__",Kingdom,sep=''),
         Phylum=paste("p__",Phylum,sep=''),
         Class=paste("c__",Class,sep=''),
         Order=if_else(is.na(Order), Order, false=paste("o__",Order,sep='')),
         Family=if_else(is.na(Family), Family, false=paste("f__",Family,sep='')),
         Genus=if_else(is.na(Genus), Genus, false=paste("g__",Genus,sep='')),
         Species=if_else(is.na(Species), Species, false=paste("s__",Species,sep='')))

#Inspect data
head(mag.taxonomy)

#Isolate and inspect tip labels. Each tip is a genomic placement
MAGtree.path <- "data/MAG.gtdbtk.bac120.classify.tree"
MAGtree <- ape::read.tree(MAGtree.path)

MAG.tips=MAGtree$tip.label
MAG.tips.keep <- MAG.tips[MAG.tips %in% mag.taxonomy$GenomeID]
tree.mags <- ape::keep.tip(MAGtree, MAG.tips.keep)


#Plot MAG tree with ggtree()
MAGtree.plot <- ggtree(tree.mags, ladderize=F) + geom_tiplab(size = 2) + xlim(0,2.6)

#Add bootstraps with the parse_bootstraps() function: 
tree.mags$node.label <- parse_bootstraps(tree.mags, method = "parse")
bs_count <- parse_bootstraps(tree.mags,method = "count")

#Finally: Add it all together for a great looking phylogeny of our MAGs:
MAGtree.plot <- MAGtree.plot + 
  geom_nodepoint(aes(subset= bs_count >= 75), fill = "cadetblue", size=2, alpha = 0.5, shape = 21) + 
  geom_treescale(x=0, y=1, width=0.1, offset = 0.5) 

mag.markingData=data.frame(GenomeID=mag.taxonomy$GenomeID, group=NA)
mag.markingData=mag.markingData %>% mutate(group=ifelse(GenomeID %in% c('metabat1500.050', 'metabat2500.019',
                                             'metabat2500.140','metabat2500.155',
                                             'concoct.124','metabat1500.117',
                                             'metabat1500.139','metabat2500.020',
                                             'metabat2500.049','metabat2500.079',
                                             'metabat2500.096','metabat2500.107',
                                             'metabat2500.120','metabat1500.113',
                                             'metabat1500.115'), 
                                    'core','other'))

MAGtree.plot=MAGtree.plot %<+% mag.markingData +
  geom_tippoint(aes(fill=group, color=group), size=3) +
  scale_fill_manual(values = c('gold', 'transparent')) +
  scale_color_manual(values = c('gold', 'transparent')) 
 

#' ************************************************************************
#' Plotting the mapping stats from samtools using MAGs and anvio contigs 
#' as template

samtoolsMAG=read.table('data/MAG_samtools_stats.txt', header=T)
samtoolsAssembly=read.table('data/anvio_samtools_stats.txt', header=T)

head(samtoolsMAG)

names(samtoolsMAG)

SamMAG.Fig=melt(samtoolsMAG[c(1,40,42)], id.vars = 'Sample', variable.name = 'portion',
     value.name = 'percent') %>%
  ggplot(aes(x=Sample, y=percent, fill=factor(portion, 
                                              levels = c('reads_unmapped_percent',
                                                         'reads_mapped_percent'))))+
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values =c('brown3', 'cadetblue'), ) +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.position = 'top')+
  labs(y='% reads', fill=NULL)


SamAssembly.Fig=melt(samtoolsAssembly[c(1,40,42)], id.vars = 'Sample', variable.name = 'portion',
                value.name = 'percent') %>%
  ggplot(aes(x=Sample, y=percent, fill=factor(portion, 
                                              levels = c('reads_unmapped_percent',
                                                         'reads_mapped_percent'))))+
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values =c('brown3', 'cadetblue'), ) +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.position = 'top')+
  labs(y='% reads', fill=NULL)


# Now plot only the MAGs without the references and include the completeness bars
# Import checkM results
MAG.checkM <- read.csv('data/MAG_checkm.csv', header = T)
head(MAG.checkM)

completness.mag.Fig <- ggplot(MAG.checkM, aes(x = fct_reorder(GeneID, -Completeness), 
                                 y = Completeness), fill='cadetblue') +
  geom_bar(stat='identity', color = "black") + #The "black" color provides the border. 
  theme_classic()+
  geom_hline(yintercept=70, color='grey70', linetype='dashed')+
  theme(legend.position = "none",
        axis.text.y= element_blank()) + #We don't need a legend for these data
  coord_flip() +
  ylab("Completeness")+
  xlab('Genome ID') #Add label 

contamination.mag.Fig <- ggplot(MAG.checkM, aes(x =fct_reorder(GeneID, -Completeness), 
                                   y = Contamination)) +
  geom_bar(stat = "identity", color = "black", fill='brown3') + #The "black" color provides the border. 
  theme_classic()+
  theme(legend.position = "none", 
        axis.text.y= element_blank(),
        axis.title.y = element_blank()) + #We don't need a legend for these data
  geom_hline(yintercept=10, color='grey70', linetype='dashed')+
  coord_flip() +
  scale_y_continuous(limits=c(0,100),oob = rescale_none)+
  ylab("Contamination")+
  xlab('Genome ID') #Add label 

grid.arrange(completness.mag.Fig, contamination.mag.Fig, widths=c(2,2))

