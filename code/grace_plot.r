#Load libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(reshape2)
library(tibble)
library(vegan)
library(BiotypeR)
library(RColorBrewer)
library(ape)
library(dendextend)
library(cowplot)
library(grid)
library(gridExtra)
library(ggdendro)
library(pairwiseAdonis)
library(stringr)

##Fig. 1A
#Plot date of collection by participant
all$date_of_supply <- as.Date(all$date_of_supply, origin = "1899-12-30")
all$date_of_sample <- as.Date(all$date_of_sample, origin = "1899-12-30")

all$max_date_collection <- with(all, ave(date_of_supply, grace_id, FUN=max))
#all$min_date_collection <- with(all, ave(date_of_supply, grace_id, FUN=min))

all$min_days <- with(all, ave(days_btwn_supply_reference, grace_id, FUN=max))

plot_date_of_supply <- ggplot(all, mapping = aes(x=date_of_supply, y=reorder(grace_id,date_of_sample)))+
  geom_point(data= na.omit(all), aes(colour=factor(drug_class)), size = 2)+
  geom_point(data=all, aes(x=date_of_sample, y=reorder(grace_id, date_of_sample), fill="Date of Sample Collection"), size=2, shape="triangle")+
  geom_segment(data= na.omit(all), aes(date_of_sample, grace_id, xend=max_date_collection, yend=grace_id))+
  #geom_hline(yintercept = 8, linetype="dotted", color = "black", size=1)+
  #geom_hline(yintercept = 16, linetype="dotted", color = "black", size=1)+
  #geom_hline(yintercept = 65, linetype="dotted", color = "black", size=1)+
  #geom_hline(yintercept = 81, linetype="dotted", color = "black", size=1)+
  scale_x_date(date_labels = "%b %y", date_breaks = "3 months")+
  theme_bw()+
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour="black"),
        legend.position = "right", legend.title = element_text(face="bold"))+
  labs(x="Date", y = "Participant", colour="Antibiotic Class", fill="Other")
plot_date_of_supply

plot_date_of_supply_colour <- plot_date_of_supply + scale_colour_manual(name="Antibiotic Class", 
                                                              values = c("Cephalosporin"="palevioletred", "Diaminopyrimidine"="lightpink", "Diaminopyrimidine and sulfonamide"="lightseagreen", 
                                                                         "Fluoroquinolone"="slategrey", "Lincosamide"="palevioletred4", "Macrolide"="darkgrey", "Methenamine Hippurate"="maroon", 
                                                                         "Nitrofuran"="peachpuff", "Nitroimidazole"="burlywood4", "Penicillin"="turquoise4", "Tetracycline"="paleturquoise"),  
                                                              breaks = c("Cephalosporin", "Diaminopyrimidine", "Diaminopyrimidine and sulfonamide", "Fluoroquinolone", "Lincosamide", 
                                                                         "Macrolide", "Methenamine Hippurate", "Nitrofuran", "Nitroimidazole", "Penicillin", "Tetracycline"),
                                                              labels = c("Cephalosporin", "Diaminopyrimidine", "Diaminopyrimidine and sulfonamide", "Fluoroquinolone", "Lincosamide", 
                                                                         "Macrolide", "Methenamine Hippurate", "Nitrofuran", "Nitroimidazole", "Penicillin", "Tetracycline"))
plot_date_of_supply_colour

##Fig. 2C
all <- read.csv("ARG_dataframe_allppts.csv", header = TRUE)

names(all)[1] <- "Description"  #Re-name first column Description

all_names <- all[,c(1,1163:1178)]
genes <- all[,1:1162]
rownames(genes) <- genes[,1]
genes <- genes[,2:1162]
sqrt_all <- sqrt(genes)

df_bray_all <- vegdist(sqrt_all, "bray")

pcoa_all <- pcoa(df_bray_all)

pcoa_variance_all <- pcoa_all$values

pcoa_all_coordinates <- data.frame(all_names,pcoa_all[["vectors"]])

gg_genes_pcoa_all <- merge(pcoa_all_coordinates,aggregate(cbind(mean.X=Axis.1,mean.Y=Axis.2)~any_tetracycline_used,pcoa_all_coordinates,mean),by="any_tetracycline_used")

pcoa_genes_all<- ggplot(data=gg_genes_pcoa_all, mapping=aes(Axis.1,Axis.2,color=any_tetracycline_used, fill = any_tetracycline_used))+ 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_text(size=11),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=8))+
  stat_ellipse(aes(Axis.1,Axis.2, color=any_tetracycline_used), geom = "polygon", alpha = 1/8, type = "t", level = 0.80)+
  theme(legend.spacing.y = unit(0.15, "cm"), legend.position = c(0.17,0.85))+
  #scale_fill_discrete(breaks=c('UTI','No UTI','Non-elderly'))+
  geom_point(size=1)+
  xlab("PCO1 (19.7% of total variation)") + ylab("PCO2 (5.9% of total variation)")+
  geom_segment(aes(mean.X, mean.Y, xend=Axis.1, yend=Axis.2))
pcoa_genes_all

##Fig. 2D
all <- read.csv("Taxa_dataframe_allppts.csv", header = TRUE)

names(all)[1] <- "Description"  #Re-name first column Description

all_names <- all[,c(1,588:603)]
taxa <- all[,1:587]
rownames(taxa) <- taxa[,1]
taxa <- taxa[,2:587]
sqrt_all <- sqrt(taxa)

df_bray_all <- vegdist(sqrt_all, "bray")

pcoa_all <- pcoa(df_bray_all)

pcoa_variance_all <- pcoa_all$values

pcoa_all_coordinates <- data.frame(all_names,pcoa_all[["vectors"]])

gg_taxa_pcoa_all <- merge(pcoa_all_coordinates,aggregate(cbind(mean.X=Axis.1,mean.Y=Axis.2)~any_tetracycline_used,pcoa_all_coordinates,mean),by="any_tetracycline_used")

pcoa_taxa_all<- ggplot(data=gg_taxa_pcoa_all, mapping=aes(Axis.1,Axis.2,color=any_tetracycline_used, fill = any_tetracycline_used))+ 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_text(size=11),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=8))+
  stat_ellipse(aes(Axis.1,Axis.2, color=any_tetracycline_used), geom = "polygon", alpha = 1/8, type = "t", level = 0.80)+
  theme(legend.spacing.y = unit(0.15, "cm"), legend.position = c(0.85,0.89))+
  #scale_fill_discrete(breaks=c('UTI','No UTI','Non-elderly'))+
  geom_point(size=1)+
  xlab("PCO1 (11.8% of total variation)") + ylab("PCO2 (6.9% of total variation)")+
  geom_segment(aes(mean.X, mean.Y, xend=Axis.1, yend=Axis.2))
  #coord_fixed()
pcoa_taxa_all
