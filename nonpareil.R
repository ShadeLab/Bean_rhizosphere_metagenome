library(Nonpareil)
library(tidyverse)

setwd('/Volumes/rs-033/ShadeLab/Stopnisek/bean_metaG/nonpareil/')
samples <- read.table('../all_samples_short.txt', sep='\t')

map <- samples %>% 
  mutate(time=if_else(grepl('pod', V1), 'pod filling', 'other'),
         time=if_else(grepl('flow', V1), 'flowering', time),
         time=if_else(grepl('senesc', V1), 'senescence', time),
         time=if_else(grepl('V2', V1), 'V2', time),
         time=if_else(grepl('V5', V1), 'V5', time),
         V1=paste(V1,".npo", sep=''),
         col=if_else(time=='V2', '#4D4D4D', ''),
         col=if_else(time=='V5', '#676767', col),
         col=if_else(time=='flowering', '#BF4942', col),
         col=if_else(time=='pod filling', '#CD6D67', col),
         col=if_else(time=='senescence', '#DA918C', col))

map_flow <- map[map$time=='senescence',]

attach(map)

np <- Nonpareil.curve.batch(V1,libnames=time, modelOnly=T, col=col, labels = time)

np_val <- c()
for(i in map$V1){
  np <- Nonpareil.curve(i)
  npVal <- np$diversity
  np_val <- rbind(np_val,npVal)
}

map_np <- cbind(map,np_val)
map_np$time <- factor(map_np$time,levels = c("V2", "V5", "flowering", "pod filling", 'senescence'))
map_np$timeNum <- if_else(map_np$time == "V2", 1, 2)

map_np$timeNum <- if_else(map_np$time == "flowering", 3, map_np$timeNum)
map_np$timeNum <- if_else(map_np$time == "pod filling", 4, map_np$timeNum)
map_np$timeNum <- if_else(map_np$time == "senescence", 5, map_np$timeNum)

nonpareilFig <-ggplot(map_np, aes(y=np_val, x=timeNum, col=time)) +
  theme_classic()+
  labs(y='Diversity Index', title='Nonpareil', x=NULL, col=NULL)+
  ylim(12,13.5)+
  geom_boxplot(outlier.colour = 'transparent') +
  geom_point(position = position_jitterdodge())+
  theme(legend.position = c(.3,.84),
        panel.background = element_rect(fill = "transparent",colour = NA))
  
Shannon_16s <- ggplot(map.alpha[map.alpha$variable == 'Shannon' & !(map.alpha$ID %in% PowersoilMRC) & map.alpha$Site == 'MRC' & map.alpha$Compartment=='rhizosphere' & map.alpha$Genotype=='Eclipse',], aes(y=value, x=factor(Timepoint), color=as.factor(Timepoint)))+
  theme_classic()+
  labs(x=NULL, y='Shannon', title='16S rRNA amplicons')+
  geom_boxplot(outlier.colour = 'transparent')+
  geom_point(position = position_jitterdodge()) +
  theme(legend.position = 'none')

ggarrange(Shannon_16s, nonpareilFig)

np <- Nonpareil.curve('MRF_senesc_3_1_S26_L007.npo')
np$diversity

detach(map)


