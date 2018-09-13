# Clear Workspace
rm(list=ls())
## Workflow Script for Sorensen_INPREP
library(ggplot2)
library(vegan)
library(outliers)
library(gplots)
library(colorRamps)
library(dplyr)
library(reshape2)
### Reading in Data Files and Manipulating Datafile
setwd("~/GitHub_Repos/ShadeLab/PAPER_Sorensen_NatMicro_2018/R")
# Mapping File
map <- read.table("Input_Files/Centralia_Collapsed_Map_forR.txt",sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)

map_MG <- map[c(1,3,4,5,6,7,10,12,14,15,16,17),]

# Carbon Dioxide data for 2015 sampling
CO2_2015 <- c(525,455,500,512,484,779,581,503,6094,17112,15760,13969,7390,723,1624,597,480,477)
CO2_MG_2015 <- CO2_2015[c(1,3,4,5,6,7,10,12,14,15,16,17)]

# 16S OTU Table
comm <- read.table("Input_Files/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)
rdp<- comm[,19]
rdp.sigs <- rdp
comm<-comm[,-19]
comm=comm[,order(colnames(comm))]

#designate a full dataset
comm.sigs=comm

#remove OTUs with an abundance = 1, across the entire dataset (singleton OTUs)
comm=comm[rowSums(comm)>1,]
sum(colSums(comm))
comm_rel <- decostand(x=comm, method="total", MARGIN=2)
colSums(comm_rel)

rdp <- rdp[rowSums(comm.sigs)>1]

# KEGG Ortholog Table
KO <- read.table("Input_Files/4-13-17abundance_ko_126107.tab.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)
colnames(KO) <- map_MG$Sample
row.names(KO) <- sub("KO:", "", row.names(KO))
KO<- KO[rowSums(KO)>0,]

##################################################################
#### Single Copy Genes abundances across temperature gradient ####
##################################################################

#Rarefy 
KO.rare <- t(rrarefy(t(KO),11855473))

#Relativized to KO count
KO.rel <- decostand(x=KO, method="total", MARGIN=2)
row.names(KO.rel) <- sub("KO:","",row.names(KO.rel))
#Single Copy KO's
library(readr)
COG_Key <- read_delim("Input_Files/He_et_al_COG_to_KEGG.txt","\t", escape_double = FALSE, trim_ws = TRUE)
COG_Key <- as.data.frame(COG_Key)
sc <- COG_Key[,3]
SCG_Rel <- KO.rel[sc,] 
SCG_Rel <- SCG_Rel[complete.cases(SCG_Rel),]
SCG_Absolute <- KO[sc,]
SCG_Absolute <- SCG_Absolute[complete.cases(SCG_Absolute),]

# For use in Odds Ratios
Average_MG_SCG <- apply(SCG_Rel, 1, mean)
#For use in Normalizing KEGG Data
Average_SCG <- apply(SCG_Absolute, 2, median)

# KOs relatvized to median Single Copy Gene abundance (This dataset will be used for further analyses)
KO.sr<- NULL
for(i in 1:nrow(KO)){
  KO.sr <- rbind(KO.sr, KO[i,]/Average_SCG)
}

KO <- as.data.frame(KO)
Odds_Ratio <- NULL
for (i in 1:nrow(SCG_Rel)){
  Odds_Ratio <- rbind(Odds_Ratio, KO.rel[row.names(SCG_Rel[i,]),]/Average_MG_SCG[i])
}

Odds_Ratio <- Odds_Ratio[complete.cases(Odds_Ratio),]


### Test for Correlations of Single Copy Gene Odds Ratios with Temperature
SCG_Correlations <- NULL
SCG_Coef <- rep(NA, nrow(Odds_Ratio))
for(i in 1:nrow(Odds_Ratio)){
  result <- cor.test(map_MG$SoilTemperature_to10cm, as.numeric(Odds_Ratio[i,]))
  SCG_Coef[i] <- lm(as.numeric(Odds_Ratio[i,])~map_MG$SoilTemperature_to10cm)$coefficients["map_MG$SoilTemperature_to10cm"]
  SCG_Correlations <- rbind(SCG_Correlations,c(unlist(result[1:4]), result$conf.int[1], result$conf.int[2]))
}
SCG_Correlations <- as.data.frame(SCG_Correlations)
SCG_Correlations$Adjusted.p.value <- p.adjust(SCG_Correlations$p.value, method="fdr")
SCG_SigCorrelations <- SCG_Correlations[SCG_Correlations$Adjusted.p.value<0.05,]

row.names(SCG_Correlations) <- row.names(Odds_Ratio)
SCG_Correlations <- as.data.frame(SCG_Correlations)
SCG_Correlations <- SCG_Correlations[order(row.names(SCG_Correlations)),]
SCG_Correlations$KEGG <- row.names(SCG_Correlations)
SCG_Correlations$RegressionCoefficients <- unlist(SCG_Coef)
SCG_Summary <- SCG_Correlations[,c(8,1,4,5,6,7,9)]
colnames(SCG_Summary) <- c("KEGG", "Test Statistic T","Pearson's r", "Lower 95% CI", "Upper 95% CI", "FDR Adjusted p-value", "Regression Coefficient" )
write.table(SCG_Summary, file = "Tables/TableS1_SingleCopyGene.txt", sep="\t", quote = FALSE)

### ggplots Odds_Rati0
MO <- cbind(Odds_Ratio, row.names(Odds_Ratio))
outl<-NULL
outa<- NULL
for (i in 1:12){
  outa <- c(outa, grubbs.test(MO[,i])[3]<0.05)
  outl <- c(outl, row.names(MO[MO[,i]==outlier(MO[,i]),])) 
}
# Remove K01519 as it is an outlier in 10/12 samples
y<-MO[MO$C12!=max(MO$C12),]
y <- melt(y, id.vars="row.names(Odds_Ratio)",variable.name="Sample", value.name="Measurement")
colnames(y) <- c("KO", "Sample", "Measurement")
# Making Data table of map and Odds Ratio values
Joined_Data <- inner_join(y, map_MG, by="Sample")
colnames(Joined_Data)[10] <- "Temperature"

z <- cor.test(Joined_Data$Temperature, Joined_Data$Measurement)

odr <- ggplot(Joined_Data, aes(x=Temperature, y=Measurement)) + 
  geom_point(size=1.5)  + 
  guides(color=FALSE) + 
  theme_bw(base_size=12) + 
  theme(text=element_text(size=8), axis.text = element_text(size=8)) +
  labs(x=expression("Temperature " ( degree~C)), y="Odds Ratio") + 
  annotate("text", x=44, y=1.4, label=paste("Pearson's r =",round(as.numeric(z[4]),3)), size=2.9) + 
  scale_y_continuous(limits=c(.5,1.5)) +
  geom_smooth(method="lm", colour="black", se=FALSE)
odr
ggsave("Figures/FigureS1A_OddsRatio.eps", odr, width=3, height=2, units="in",device="eps" )

########################################
#### Look at Eukaryotic DNA content ####
########################################
Eukaryotes <- read.table("Input_Files/Eukaryotic_Ribosome_Simple.txt", sep="\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)

KO.rare.eu <- KO.rare[Eukaryotes[Eukaryotes$Taxonomy=="Eukaryotes",1],]
Eukarytoes_Correlation <- data.frame(KO=row.names(KO.rare.eu), t_value=rep(NA,nrow(KO.rare.eu)), Pearsons_R=rep(NA, nrow(KO.rare.eu)), p_value= rep(NA, nrow(KO.rare.eu)))
for(i in 1:nrow(KO.rare.eu)){
  x <- cor.test(KO.rare.eu[i,], map_MG$SoilTemperature_to10cm)
  Eukarytoes_Correlation$Pearsons_R[i] <- x$estimate
  Eukarytoes_Correlation$p_value[i] <- x$p.value
  Eukarytoes_Correlation$t_value[i] <- x$statistic
}

write.table(x=Eukarytoes_Correlation, file = "Tables/EukaryoticSCG_Correlations.txt", sep="\t", quote = FALSE, row.names = FALSE)

#############################################################
#### Calculate Each Phylum's Average Genome Size in IMG  ####
#############################################################
JGI_Archaea <- read.table("Input_Files/JGI_Archaea_06192017.txt", sep="\t", header=TRUE, row.names=NULL, stringsAsFactors = FALSE, fill=TRUE, quote="")

Permanent_Archaea <- JGI_Archaea[grepl("Permanent", JGI_Archaea$Status),]
Finished_Archaea <- JGI_Archaea[grepl("Finished",JGI_Archaea$Status),]


JGI_Bacteria <- read.table("Input_Files/JGI_Bacteria_06192017.txt", sep="\t", header=TRUE, row.names=NULL, stringsAsFactors = FALSE, fill=TRUE, quote="")

Permanent_Bacteria <- JGI_Bacteria[grepl("Permanent", JGI_Bacteria$Status),]
Finished_Bacteria <- JGI_Bacteria[grepl("Finished", JGI_Bacteria$Status),]

PF_Archaea <- rbind(Permanent_Archaea, Finished_Archaea)
PF_Bacteria<- rbind(Permanent_Bacteria,Finished_Bacteria)

PF_All <- rbind(PF_Archaea, PF_Bacteria)

#Formatting Taxonomy from 16S data
rdp <- rdp[rowSums(comm.sigs)>1]
rdp <- gsub("\\[|\\]","", rdp)

Taxonomy <- NULL
for(i in 1:length(rdp)){
  Taxonomy <- rbind(Taxonomy, unlist(strsplit(as.character(rdp[i]), ";")))
}

Phylogeny <- Taxonomy[,2]
Phyla <- Taxonomy[,2]
Phyla<- gsub(" p__", "", Phyla)
Phyla <- gsub("k__", "", Phyla)
Phyla <- unique(Phyla)
Phyla[which(!Phyla%in%PF_All$Phylum)]

Phylogeny <- Taxonomy[,2]

# Detected taxa with no sequenced genomes in IMG. Call them Bacteria
Not_Present <- c("AC1", "AD3", "FCPU426", "GN02", "GN04", "GOUTA4", "LD1", "MVP-21", "MVS-104", "PAUC34f", "SAR406", "SBR1093", "SC4", "WS2","WS3", "WS4", "WS5", "WWE1")
for(i in 1:length(Not_Present)){
  Phylogeny<-gsub(Not_Present[i], "Bacteria", Phylogeny )
}
# Key for converting Greengenes taxonomy to JGI/IMG taxonomy
GG_to_JGI <- read.table("Input_Files/GG_to_JGI.txt", header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE)

for(i in 1:nrow(GG_to_JGI)){
  Phylogeny <- gsub(GG_to_JGI[i,1], GG_to_JGI[i,2], Phylogeny)
}

Phylogeny<- gsub(" p__", "", Phylogeny)
Phylogeny <- gsub("k__", "", Phylogeny)
v <- unique(Phylogeny)
v<- v[-42]
Fixed_Phylum <- unique(Phylogeny)
Fixed_Phylum[22] <- "Candidatus Parvarchaeota"
Fixed_Phylum <- Fixed_Phylum[complete.cases(Fixed_Phylum)]
Fixed_Phylum <- Fixed_Phylum[-42]  


# For each phylum calculate
# Total,Average, IQR(Q3-Q1), Q1, Q3, Lower inner fence(Q1-1.5*IQR), Upper inner fence(Q3+1.5*IQR), Outliers, New Mean, SD  

Necessary_JGI <- NULL
for (i in 1:length(Fixed_Phylum)){
  x <- PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14]
  Total <- as.numeric(length(x))
  avg_bad <- mean(x) 
  Q1<- as.numeric(quantile(x)[2])
  Q3 <- as.numeric(quantile(x)[4])
  IQR <- Q3-Q1
  LIF <- Q1-(1.5*IQR)
  UIF <- Q3 + (1.5*IQR)
  O <- Total-sum(1*(LIF<x & x<UIF))# Number of outlier genomes for that phylum
  x_good <- x[LIF<x & x<UIF]
  avg_good <- mean(x_good)
  sd_good <- sd(x_good)
  Necessary_JGI <- rbind(Necessary_JGI, c(Total, avg_bad, Q1, Q3, IQR, LIF, UIF, O, avg_good, sd_good, sd_good/avg_good))
}
row.names(Necessary_JGI)<- Fixed_Phylum

Necessary_JGI2 <- Necessary_JGI
colnames(Necessary_JGI2) <- c("Total_Genomes", "Average_Total", "Q1", "Q3","IQR","LIF", "UIF", "Outliers", "Average_Clean", "SD_Clean", "CV")

Necessary_JGI2 <- as.data.frame(Necessary_JGI2)

TableSX <- data.frame(Phylum = row.names(Necessary_JGI2), Genomes=Necessary_JGI2$Total_Genomes - Necessary_JGI2$Outliers, AverageGenomeSize = Necessary_JGI2$Average_Clean)

write.table(file="Tables/TableSX_PhylumGenomeSize.txt", TableSX, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

plot_data<- NULL
for (i in 1:length(Fixed_Phylum)){
  x <- cbind(rep(Fixed_Phylum[i], length(PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14])), PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14])
  plot_data <- rbind(plot_data, x)             
}
plot_data <- as.data.frame(plot_data)



x <- PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14]

###########################################################################
#### Average Genome Size based on Phylum Level Abundance from 16S data ####
###########################################################################

#creating Phylum level table for Centralia Data
comm_12<-comm[,c(1,3,4,5,6,7,10,12,14,15,16,17)]
Phyla_Comm<- NULL
for(i in 1:length(Fixed_Phylum)){
  x <- comm_12[grepl(Fixed_Phylum[i],Phylogeny),]
  y<- as.numeric(colSums(x))
  Phyla_Comm <- rbind(Phyla_Comm, y)
}
row.names(Phyla_Comm) <- Fixed_Phylum
Phyla_Comm_Rel <- decostand(Phyla_Comm, method= "total", MARGIN=2 )
colSums(Phyla_Comm_Rel)




sum(1*(PF_Bacteria$Genome.Size.....assembled<100000))
nrow(PF_Bacteria)
Necessary_JGI["Bacteria", 9] <- mean(PF_Bacteria$Genome.Size.....assembled)
Necessary_JGI["Unclassified", 9] <- mean(PF_All$Genome.Size.....assembled)
Necessary_JGI["OP11", 9] <- mean(PF_Bacteria$Genome.Size.....assembled)
Necessary_JGI["FBP", 9] <- mean(PF_Bacteria$Genome.Size.....assembled)
Necessary_JGI["Archaea", 9] <- mean(PF_Archaea$Genome.Size.....assembled)
Necessary_JGI["candidate division TM6", 9] <- Necessary_JGI["candidate division TM6",2]
Necessary_JGI["WS1", 9] <- Necessary_JGI["WS1",2]
Necessary_JGI["candidate division SR1", 9] <- Necessary_JGI["candidate division SR1",2]


Output2 <- NULL
for(i in 1:ncol(Phyla_Comm_Rel)){
  x <- Phyla_Comm_Rel[,i]*Necessary_JGI[,9]
  z <- sum(x)
  Output2 <- c(Output2, z)
}

plot(map_MG$SoilTemperature_to10cm, Output2)
cor.test(map_MG$SoilTemperature_to10cm, Output2)
plot_data_2 <- cbind(Output2, map_MG$SoilTemperature_to10cm)
plot_data_2 <- as.data.frame(plot_data_2)
colnames(plot_data_2) <- c("Estimate", "Temperature")

#### Plot 16S rRNA Estimate
z <- cor.test(plot_data_2$Temperature, plot_data_2$Estimate)
rRNA_Size_Estimate <- ggplot(plot_data_2, aes(x=Temperature, y =Estimate )) + 
  geom_point(size=1.5) + 
  theme_bw(base_size=8) + 
  theme(text=element_text(size=8), axis.text = element_text(size=8)) +
  labs(x=expression("Temperature " ( degree~C)), y="Average Genome Size (bp)") +
  annotate("text", x=45, y=3600000, label=paste("Pearson's r =",round(as.numeric(z[4]),3)), size=2.9) +
  scale_y_continuous(labels = function(y) format(y, scientific = TRUE), limits=c(3000000, 3800000))+
  geom_smooth(method="lm", colour="black", se=FALSE)
rRNA_Size_Estimate
ggsave("Figures/FigureS1B_rRNA_Estimate.eps", rRNA_Size_Estimate, width=3, height=2, units="in",device="eps" )

######################################################
#### Average Genome Sizes based on Microbe Census ####
######################################################

# Microbe Census Estimates
M_Census <- c(6705016, 6046631, 6174721, 6213575, 6088619, 5810701, 4136147, 5295248, 4506566, 4881807, 5421668, 6144962)
plot_data_3 <- cbind(M_Census, map_MG$SoilTemperature_to10cm)
plot_data_3 <- as.data.frame(plot_data_3)
colnames(plot_data_3) <- c("Estimate", "Temperature")

z <- cor.test(plot_data_3$Temperature, plot_data_3$Estimate)
MCensus_Size_Estimate <- ggplot(plot_data_3, aes(x=Temperature, y =Estimate )) + 
  geom_point(size=1.5) + 
  theme_bw(base_size=8) + 
  theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x=expression("Temperature " ( degree~C)), y="Average Genome Size (bp)")+ 
  annotate("text", x=45, y=6000000, label=paste("Pearson's r =",round(as.numeric(z[4]),3)), size=2.9) + 
  scale_y_continuous(labels = function(y) format(y, scientific = TRUE),limits=c(4000000, 7000000)) +
  geom_smooth(method="lm", color="black", se=FALSE)

ggsave("Figures/Figure1A_MCensus.eps", MCensus_Size_Estimate, width=3, height=2, device="eps", units="in")

### Multiplot code taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##############################
#### KEGG Module Analysis ####
##############################
library(stringr)
Modules <- read.table("Input_Files/kmodlist47982_23-nov-2016.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

# Find those KEGG Modules which have other KEGG Modules as part of their definitions
MO_Per_M <- rep(0,nrow(Modules))
for (i in 1:nrow(Modules)){
  MO_Per_M[i] <- str_count(Modules$Definition[i],"M")
}

Modules_w_Modules <- Modules[MO_Per_M>0,]

Dictionary=NULL
Dictionary <- vector(mode="list", length=nrow(Modules))
names(Dictionary) <- row.names(Modules)

Dictionary.rel=NULL
Dictionary.rel <- vector(mode="list", length=nrow(Modules))
names(Dictionary.rel) <- row.names(Modules)

### For Some Module definitions, they have other modules in them, this is formating those definitions so that there are only K0 and no M0's in the Definition. As a side note these are all "Signature" module types
#Acetogen
Modules["M00618",3] <- paste(c(Modules["M00377",3],Modules["M00579",3]), collapse=" ")
#Anoxygenic photosynthesis in green nonsulfur bacteria
Modules["M00613",3] <- paste(c(Modules["M00597",3],Modules["M00376",3]), collapse=" ")
#Anoxygenic photosynthesis in green sulfur bacteria
Modules["M00614",3] <- paste(c(Modules["M00598",3],Modules["M00173",3]), collapse=" ")
#Anoxygenic photosynthesis in purple bacteria
Modules["M00612",3] <- paste(c(Modules["M00597",3],Modules["M00165",3]), collapse=" ")
#Methanogen
Modules["M00617",3] <- paste(c(Modules["M00567",3],Modules["M00357",3],Modules["M00356",3],Modules["M00563",3]), collapse=" ")
#Nitrate assimilations
Modules["M00615",3] <- paste(c(Modules["M00438",3],Modules["M00531",3]), collapse=" ")
#Oxygenic photosynthesis in plants and cyanobacteria 
Modules["M00611",3] <- paste(c(Modules["M00161",3],Modules["M00163",3],Modules["M00165",3]), collapse=" ")
#Sulfate-sulfur assimilation
Modules["M00616",3] <- paste(c(Modules["M00185",3],Modules["M00176",3]), collapse=" ")

#Produces a list, each item in this list is a dataframe with the KO's in our dataset 
# Relativized to Average Single Copy Gene count
for (y in 1:nrow(KO.sr)){
  KO_M <- grep(row.names(KO.sr)[y], Modules$Definition)
  if (length(KO_M)>0){
    for (x in 1:length(KO_M)){
      Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.sr[y,])
    }
  }
}

### Number of modules not represented by dataset
count <- 0
for (y in 1:length(Dictionary)){
  if (data.class(Dictionary[[y]])=="NULL"){
    count <- count + 1
  }
}
count

#Counts the number of KOs in each module
library(stringr)
KO_Per_M <- rep(0, nrow(Modules))
for (i in 1:nrow(Modules)){
  KO_Per_M[i] <- str_count(Modules$Definition[i],"K")
}
# Counts number of KOs from each module that are present in our dataset
KO_Per_M_Data <- rep(0, nrow(Modules))
for (i in 1:nrow(Modules)){
  if(data.class(Dictionary[[i]])!="NULL"){
    KO_Per_M_Data[[i]] <- nrow(Dictionary[[i]])
  }
}

Missing_KOs <- KO_Per_M - KO_Per_M_Data

KO_Per_M_Present <- KO_Per_M[KO_Per_M_Data>0]
KO_Per_M_Data_Present <- KO_Per_M_Data[KO_Per_M_Data>0]
### Number of Modules that are less than 50% complete 
sum(1*((KO_Per_M_Data_Present/KO_Per_M_Present)<0.5))

length(Missing_KOs[Missing_KOs>0])

KO_Per_M_Data/KO_Per_M


### Getting Rid of Modules that are not present at all in our data. 
Dictionary_Subset <- Dictionary[lapply(Dictionary, data.class) == "data.frame"]

Modules_Subset <- Modules[lapply(Dictionary,data.class) == "data.frame",]
KO_Per_M_Subset <- KO_Per_M[lapply(Dictionary, data.class) == "data.frame"]
KO_Per_M_Data_Subset <- KO_Per_M_Data[lapply(Dictionary, data.class) == "data.frame"]

Dictionary.rel_Subset <- Dictionary.rel[lapply(Dictionary.rel, data.class)=="data.frame"]


#Fraction of KOs from each module present in our dataset
Module_Completeness <- KO_Per_M_Data_Subset/KO_Per_M_Subset


#Calculating correlation coefficient and T.test results for all KOs in all modules  
Summary_Output <- NULL
avg_mod <- NULL
mid_mod <- NULL
max_mod <- NULL
min_mod <- NULL
stdev_mod <- NULL
for (m in 1:length(Dictionary_Subset)){
  temp_data <- Dictionary_Subset[[m]]
  K_out <- NULL
  agg <- NULL
  ### Calculates the AVG, Median, Min, Max, and STDEV for the KOs of a module in each given site.
  for (s in 1:ncol(temp_data)){
    avg <- mean(temp_data[,s])
    mid <- median(temp_data[,s])
    minimum <- min(temp_data[,s])
    maximum <- max(temp_data[,s])
    stdev <- sd(temp_data[,s])
    agg <- cbind(agg, c(avg, mid, minimum, maximum, stdev))
  }
  avg_mod <- rbind(avg_mod, agg[1,])
  mid_mod <- rbind(mid_mod, agg[2,])
  min_mod <- rbind(min_mod, agg[3,])
  max_mod <- rbind(max_mod, agg[4,])
  stdev_mod <- rbind(stdev_mod, agg[5,])
}

colnames(mid_mod) <- map_MG$Sample
row.names(mid_mod) <- names(Dictionary_Subset)

Module_Completeness_Final <- Module_Completeness[names(Dictionary_Subset)%in%row.names(mid_mod)]


KO_Per_M_Data_Present <- KO_Per_M_Data_Present[rowSums(mid_mod)!=0]
KO_Per_M_Present <- KO_Per_M_Present[rowSums(mid_mod)!=0]
mid_mod <- mid_mod[rowSums(mid_mod)!=0,]

Module_Completeness_Final <- Module_Completeness[names(Dictionary_Subset)%in%row.names(mid_mod)]

mid_mod_PA<- 1*(mid_mod>0) 
#Percent of modules present in every site
sum(1*(rowSums(mid_mod_PA)==12))/nrow(mid_mod_PA)

### *** Module Temperature Correlations, T-Tests, and Heatmaps *** 
median_Module_TempCorrelations <- apply(mid_mod, 1, function(x) cor.test(x, map_MG$SoilTemperature_to10cm))

Reg_Coef<- rep(NA,nrow(mid_mod))
for(i in 1:nrow(mid_mod)){
  Reg_Coef[i] <- lm(mid_mod[i,]~map_MG$SoilTemperature_to10cm)$coefficients["map_MG$SoilTemperature_to10cm"]
}
### 229 Modules significantly correlated with temperature
Med_Mod_TempCor <- NULL
for (i in 1:length(median_Module_TempCorrelations)){
  Med_Mod_TempCor <- rbind (Med_Mod_TempCor, c(median_Module_TempCorrelations[[i]][4], median_Module_TempCorrelations[[i]]$conf.int[1], median_Module_TempCorrelations[[i]]$conf.int[2], median_Module_TempCorrelations[[i]][2],median_Module_TempCorrelations[[i]][3]))
}
Med_Mod_TempCor <- as.data.frame(Med_Mod_TempCor)
Med_Mod_TempCor$KEGG <- row.names(mid_mod)
Med_Mod_TempCor[,1] <- unlist(Med_Mod_TempCor[,1])
Med_Mod_TempCor[,2] <- unlist(Med_Mod_TempCor[,2])
Med_Mod_TempCor[,3] <- unlist(Med_Mod_TempCor[,3])
Med_Mod_TempCor[,4] <- unlist(Med_Mod_TempCor[,4])
Med_Mod_TempCor[,5] <- unlist(Med_Mod_TempCor[,5])
Med_Mod_TempCor$Completeness <- KO_Per_M_Data_Present/KO_Per_M_Present
Med_Mod_TempCor$Adjusted.p.value <- p.adjust(Med_Mod_TempCor$p.value, "fdr")
row.names(Med_Mod_TempCor)<- Med_Mod_TempCor$KEGG
Med_Mod_TempCor$RegressionCoefficient <- Reg_Coef


Modules_Subset <- Modules_Subset[row.names(Med_Mod_TempCor),]

Modules_Subset$Completeness <- Module_Completeness_Final

Complete_Modules<- Modules_Subset[]

Med_Mod_TempCor <- Med_Mod_TempCor[Med_Mod_TempCor$Completeness>=0.5,]

Sig_Temp_Cor <- Med_Mod_TempCor[Med_Mod_TempCor$Adjusted.p.value<0.05,]

colnames(Sig_Temp_Cor)[2] <- "Lower 95% CI"
colnames(Sig_Temp_Cor)[3] <- "Upper 95% CI"

Modules_Subset_Complete <- Modules_Subset[Module_Completeness_Final>=0.5,]
#library(plyr)
Combined_Sig_Module_Results <- Sig_Temp_Cor
Combined_Sig_Modules <- as.data.frame(mid_mod)[Combined_Sig_Module_Results$KEGG,]

### Create Supplemental Table of Summaries
SignificantModules_Summary <- Combined_Sig_Module_Results

SignificantModules_Summary$ModuleDescription <- Modules[row.names(Combined_Sig_Modules),1]

row.names(SignificantModules_Summary) <- SignificantModules_Summary$KEGG 

SignificantModules_Summary<-SignificantModules_Summary[,c(6,10,7,1,2,3,9,8)]
colnames(SignificantModules_Summary) <- c("Module","Module Description","Completeness","Pearson's r", "Lower 95% CI", "Upper 95% CI","Regression Coefficient", "FDR adjusted p-value")
write.table(x = SignificantModules_Summary, file="Tables/TableS3_Median_SCG.txt", sep="\t", quote=FALSE)


### Heat map of Modules that are either significant based in T-Test or Temperature Correlation
library(colorRamps)
library(gplots)
library(vegan)
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

# Two-component Regulatory System Results
TCRS <- SignificantModules_Summary[grep("two-component", SignificantModules_Summary$`Module Description`),]
TCRS_M <- Combined_Sig_Modules[TCRS[,1],]
med_TCRS <- apply(TCRS_M, 2, median)

# Drug resistance Results
All_Drug <- read.table("Input_Files/All_Drug_Modules.txt", sep="\t", row.names=NULL, header=TRUE, stringsAsFactors = FALSE)

All_Drug_Stats <- Med_Mod_TempCor[All_Drug$Module.ID,]
All_Drug_Stats <- All_Drug_Stats[complete.cases(All_Drug_Stats),]


Combined_Sig_Modules.zs <- decostand(as.matrix(Combined_Sig_Modules), method="standardize", MARGIN=1)
row.names(Combined_Sig_Modules.zs)

Pos_Temp_Cor_Modules <- Combined_Sig_Modules[Combined_Sig_Module_Results$estimate>0,]
#Positively Correlated Modules z-scored
PCM.zs <- decostand(as.matrix(Pos_Temp_Cor_Modules), method="standardize", MARGIN=1)

PCM.zs <- PCM.zs[,order(map_MG$SoilTemperature_to10cm)]

# Negative Temp Cor Modules
Neg_Temp_Cor_Modules <- Combined_Sig_Modules[Combined_Sig_Module_Results$estimate<0,]
# Z-score Negative Correlated Modules
NCM.zs <- decostand(as.matrix(Neg_Temp_Cor_Modules), method="standardize", MARGIN = 1)
NCM.zs <- NCM.zs[,order(map_MG$SoilTemperature_to10cm)]


# Plot and Save PCM heatmap Key
setEPS()
postscript("Figures/Figure4A_PositiveCorrelatedHeatmap.eps", width = 3.4, height=9, paper="special")
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(PCM.zs, col=hc(100), key=TRUE, symkey=TRUE, trace="none",density.info = "none", colsep=c(1:12),rowsep=c(1:nrow(PCM.zs)), sepcolor="black", Colv=FALSE, sepwidth=c(0.01,0.00001),labRow = row.names(PCM.zs), dendrogram="row",cexRow = 1, srtCol=90, lmat= rbind(c(3,4),c(2,1)), lhei=c(1,4))
dev.off()

NCM_Row_Colors<- rep(NA, nrow(NCM.zs))
NCM_Row_Colors[row.names(NCM.zs)%in%row.names(TCRS)] <- "blue"
NCM_Row_Colors[row.names(NCM.zs)%in%All_Drug_Stats$KEGG] <- "grey"
NCM_Row_Colors[row.names(NCM.zs)%in%All_Drug_Stats$KEGG&row.names(NCM.zs)%in%row.names(TCRS)] <- "black"


# Plot and Save NCM heatmap
setEPS()
postscript("Figures/Figure4B_NegativeCorrelatedHeatmap.eps", width = 3.4, height=9, paper="special")
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(NCM.zs, col=hc(100),key=FALSE, symkey=TRUE, trace="none", colsep=c(1:12),density.info="none", sepcolor="black", Colv=FALSE, sepwidth=c(0.01,0.0000001), dendrogram="row",cexRow = 1, labRow=FALSE, srtCol=90,lhei=c(1,100), RowSideColors = NCM_Row_Colors)
dev.off()
par(mfrow=c(1,1))

# Plot Temperature Circles for Heatmaps
map_MG$Same <- rep(1,nrow(map_MG))
GnYlOrRd=colorRampPalette(colors=c("green", "yellow", "orange","red"), bias=2)
Temperature_Points <- ggplot(map_MG, aes(x=map_MG$Sample, y=Same)) + geom_point(aes(y=as.numeric(Same), x=factor(map_MG$Sample, levels=map_MG$Sample[order(map_MG$SoilTemperature_to10cm)]), color=as.numeric(SoilTemperature_to10cm)), size=5) + scale_color_gradientn(colours=GnYlOrRd(5), guide="colorbar", guide_legend(title="Temperature")) +theme_bw()

ggsave("Figures/Figure4_TemperatureCircles.eps",width = 4, units = "in", height = 4, dpi=300 )

#######################################
#### Calculating average cell size ####
#######################################
Cell_Size_Names <- c("Cen01", "Cen03", "Cen04", "Cen05", "Cen06", "Cen07", "Cen10", "Cen12", "Cen14", "Cen15", "Cen16", "Cen17")
Cell_Size_Summary <- matrix(nrow=12, ncol=7)
m <- NULL
n<-NULL
q <- NULL
for (i in 1:length(Cell_Size_Names)){
  x <- read.table(paste("Input_Files/CleanImages/", Cell_Size_Names[i],".csv", sep=""), sep=",",header = TRUE, row.names = 1 )
  Cell_Size_Summary[i,] <- c(mean(x[,1]), sd(x[,1]), mean(x[,7]), sd(x[,7]), mean(x[,8]), sd(x[,8]), nrow(x))
  m <- rbind(m, x[1:nrow(x),])
  n <- c(n,rep(Cell_Size_Names[i], nrow(x)))
}

Cell_Size_Summary
row.names(Cell_Size_Summary) <- Cell_Size_Names
colnames(Cell_Size_Summary) <- c("Average_Area", "SD_Area", "Average_Length", "SD_Length","Average_Minor", "SD_Minor", "Number_Cells")
Cell_Size_Summary <- as.data.frame(Cell_Size_Summary)

write.table(file="Tables/TableS4_CellSizeSummary.txt",x=Cell_Size_Summary, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
Cell_Size_Summary$Temperature <- map_MG$SoilTemperature_to10cm

Cell_Size_Cor <- cor.test(Cell_Size_Summary$Average_Length,Cell_Size_Summary$Temperature)
Cell_Length_Plot <- ggplot(Cell_Size_Summary, aes(x=Temperature, y=Average_Length)) + 
  geom_point(size=1.5)  + 
  guides(color=FALSE) + 
  theme_bw(base_size=12) + 
  theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x=expression("Temperature " ( degree~C)), y=expression(paste("Average Length (", mu ,"m)"))) +
  scale_y_continuous(limits=c(.7, 1))+
  annotate("text", x=45, y=0.91, label=paste("Pearson's r =",round(as.numeric(Cell_Size_Cor[4]),3)), size=2.9)+
  geom_smooth(method="lm", color="black", se=FALSE)

Cell_Length_Plot
ggsave("Figures/Figure1B_CellLength.eps", width=3, height=2, units="in", device ="eps" )


cor.test(Cell_Size_Summary$Average_Length, M_Census)

Cell_Size_Summary$Genome_Size <- M_Census
Length_v_Size_Cor <- cor.test(Cell_Size_Summary$Average_Length, M_Census)

Length_v_Genome_Plot <- ggplot(Cell_Size_Summary, aes(x=Genome_Size, y=Average_Length)) + 
  geom_point(size=1.5)  + 
  guides(color=FALSE) + 
  theme_bw(base_size=12) + 
  theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x="Average Genome Size (bp)", y=expression(paste("Average Length (", mu ,"m)"))) +
  scale_y_continuous(limits=c(.7, 1))+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE),limits=c(4000000, 7000000))+
  annotate("text", x=4750000, y=0.91, label=paste("Pearson's r =",round(as.numeric(Length_v_Size_Cor[4]),3)), size=2.9)+
  geom_smooth(method="lm", color="black", se=FALSE)

ggsave("Figures/Figure1C_Length_v_Genome.eps", Length_v_Genome_Plot, width=3, height=2, device="eps", units="in")


colnames(Cell_Size_Summary) <- c("Average_Area", "SD_Area", "Average_Length", "SD_Length","Average_Minor", "SD_Minor", "Number_Cells")
row.names(Cell_Size_Summary) <- Cell_Size_Names

##########################
#### Average MAG Size ####
##########################
b <- read.table("Input_Files/CheckM_VerySpecific_Results.txt", stringsAsFactors = FALSE)
b <- b[,-3]
colnames(b)<- c("BinID", "MarkerLineage", "#Genomes", "#Markers", "#MarkerSets", "0","1", "2", "3", "4", "5+", "Completeness", "Contamination", "StrainHeterogeneity")

# 90% or Greater Completeness
Complete <- b[b$Completeness>=90,]
Complete_Clean <- Complete[Complete$Contamination<5,]
write.table(Complete_Clean$BinID, file="~/90Complete_5Contam.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

#write.table(Complete_Clean[,1:13], file="Tables/TableS6_MAG_Details.txt", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)

# Read in Coverage Files for Each MAG
Files <- system(command = "ls Input_Files/MAG_Coverage/Coverage*", intern=TRUE)



MAGS <- rep("A", length(Files))
for (i in 1:length(Files)){
  x <- unlist(strsplit(Files[i], "/"))[3]
  MAGS[i] <- gsub(x,pattern="Coverage.", replacement= "")
}


MAGS_Coverage <- matrix(nrow=length(Files), ncol=12)
MAGS_Length <- rep(1, length(Files))
for(i in 1:length(Files)){
  x <- read.table(Files[i], sep="\t", header=TRUE, row.names=NULL, stringsAsFactors = FALSE)
  MAGS_Length[i] <- sum(x$contigLen)
  MAGS_Coverage[i,] <- apply(x[,c(4,6,8,10,12,14,16,18,20,22,24,26)], MARGIN=2, mean)
}
MAGS_Coverage
row.names(MAGS_Coverage) <- MAGS
MAG_GenomeSize_Estimate <- MAGS_Length*(100/Complete_Clean$Completeness)
MAG_GenomeSize_Estimate <- data.frame(MAG=MAGS, Genome_Length=MAG_GenomeSize_Estimate)

hist(MAG_GenomeSize_Estimate$Genome_Length)


Site_Average_Genome_Size_MAG <- rep(0,12)
for (i in 1:ncol(MAGS_Coverage[,1:12])){
  Site_Average_Genome_Size_MAG[i] <- mean(MAG_GenomeSize_Estimate$Genome_Length[MAGS_Coverage[,i]>0])
}

MAG_Plot_Data <- data.frame(Size=Site_Average_Genome_Size_MAG, Temperature=map_MG$SoilTemperature_to10cm)
MAG_Cor <- cor.test(map_MG$SoilTemperature_to10cm, Site_Average_Genome_Size_MAG)


ggplot(MAG_Plot_Data, aes(x=Temperature, y=Size)) + 
  geom_point(size=1.5)  + 
  guides(color=FALSE) + 
  theme_bw(base_size=12) + 
  theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x=expression("Temperature " ( degree~C)), y="Average MAG Size (bp)")+
  scale_y_continuous(labels = function(y) format(y, scientific = TRUE), limits=c(4200000, 4600000))+
  annotate("text", x=45, y=4550000, label=paste("Pearson's r =",round(as.numeric(MAG_Cor[4]),3)), size=2.9)+
  geom_smooth(method="lm", color="black", se=FALSE)

ggsave("Figures/FigureS1C_AVG_MAG_Size.eps", width=3, height=2, units="in", device ="eps" )

#load packages
library(reshape2)
library(tidyverse)
library(psych)


##View trend of average genome size 
#read in microbe census data
Fig2Data <- read.delim("Input_Files/AGS_cent_references.txt")

#read in site information
meta <- read.delim(file = "Input_Files/sample_map.txt", header = TRUE)

#join datasets to annotate microbe census data
Fig2Data.annotated <- meta %>%
  dplyr::rename(Project = Site, Site = Sample) %>%
  mutate(Site = gsub("cen", "Cen", Site)) %>%
  right_join(Fig2Data, by="Site") %>%
  mutate(Site = gsub("_", " ", Site),
         Project = gsub("_", " ", Project))

#remove samples that were eliminated
Fig2Data.annotated <- Fig2Data.annotated[!is.na(Fig2Data.annotated$Biome),]

#plot data
(AGS <- ggplot(Fig2Data.annotated, aes(x=reorder(Site, AGS), y=AGS/1000000)) +
    geom_col(aes(fill = Project), color = "black") +
    ylab("Average genome size (Mbp)") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = c("grey65", "red", "yellow1", "green3","#E495A5", "#CBA56F", "#B2AE64", "#92B66E","#6DBC86", "#46BEA3", "#3ABCBE", "#5EB6D4", "#8CADE0", "#B5A1E1","#D198D6", "grey90")) +
    xlab("Microbial Samples") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

#save plot
ggsave(AGS, filename = "Figures/Figure2_AGS_Barplots.eps", height = 6, width = 7, units = "in")
AGS

################################################
#### HOBO Data for temperature supplemental ####
################################################
Cen08_2016 <- read.table("Input_Files/HOBO_Data/Cen08_12102016_21102017_0.txt", sep="\t", header=TRUE, row.names=1)

Cen08_2016$Date <- as.Date(Cen08_2016$Date, format="%m-%d-%y")

Cen23_2016 <- read.table("Input_Files/HOBO_Data/Cen23_13102016_21102017.txt", sep="\t", header=TRUE, row.names=1)

Cen23_2016$Date <- as.Date(Cen23_2016$Date, format="%m-%d-%y")

Cen19_2016 <- read.table("Input_Files/HOBO_Data/Cen19_12102016_19102017.txt", sep="\t", header=TRUE, row.names=1)

Cen19_2016$Date <- as.Date(Cen19_2016$Date, format="%m-%d-%y")

Cen14_2016 <- read.table("Input_Files/HOBO_Data/Cen14_11102016_21102017.txt", sep="\t", header=TRUE, row.names=1)
hist(Cen14_2016$Temp...C..LGR.S.N..10809340..SEN.S.N..10809340.)
# Let's clean up this temp data, look like most temp readings are between 0 and 50C, lets remove any <-10Cor >75C
Cen14_2016 <- Cen14_2016[Cen14_2016$Temp...C..LGR.S.N..10809340..SEN.S.N..10809340.>0,]
Cen14_2016 <- Cen14_2016[Cen14_2016$Temp...C..LGR.S.N..10809340..SEN.S.N..10809340.<50,]

Cen14_2016$Date <- as.Date(Cen14_2016$Date, format="%m-%d-%y")
range(Cen14_2016$Date)


Cen15_2016 <- read.table("Input_Files/HOBO_Data/Cen15_11102016.txt", sep="\t", header=TRUE, row.names=1)
Cen15_2016$Date <- as.Date(Cen15_2016$Date, format="%m-%d-%y")
range(Cen15_2016$Date)


Dates <- unique(Cen14_2016$Date)
v <- rep(NA, length(Dates))
w <- rep(NA, length(Dates))
x <- rep(NA, length(Dates))
y <- rep(NA, length(Dates))
z <- rep(NA, length(Dates))
for(i in 1:length(Dates)){
  v[i] <- mean(Cen08_2016[Cen08_2016$Date==Dates[i], 3])
  w[i] <- mean(Cen14_2016[Cen14_2016$Date==Dates[i], 3]) 
  x[i] <- mean(Cen15_2016[Cen15_2016$Date==Dates[i], 3]) 
  y[i] <- mean(Cen19_2016[Cen19_2016$Date==Dates[i], 3])
  z[i] <- mean(Cen23_2016[Cen23_2016$Date==Dates[i], 3])
}

Plot_Data <- data.frame(Date=Dates, Cen08=v, Cen14=w, Cen15=x, Cen19=y, Cen23=z )

Plot_Data_m <- melt(Plot_Data, id.vars = "Date", measure.vars = c("Cen08","Cen14","Cen15","Cen19", "Cen23"))
colnames(Plot_Data_m)[2] <- "Site"

Plot_Data_m$Classification <- c(rep("Reference", length(Dates)), rep("FireAffected", 3*length(Dates)), rep("Reference", length(Dates)))


ggplot(Plot_Data_m, aes(x=Date, y=value)) +
  geom_point(aes(shape=Classification, color=Site )) +
  labs(y=expression(" AVerage Daily Temperature " ( degree~C)))

ggsave("Figures/FigureS2_HOBO_2016_2017.eps", device="eps", width=6, height=4, units="in")


###############################################################
#### Contextual Data correlations with Average Genome Size ####
###############################################################

Contextual_data <- map_MG[,c(8,13:24)]
Contextual_data_summary <- data.frame(Pearsons_R=rep(NA,12), Test_Statistic=rep(NA,12), p_value=rep(NA,12))
for(i in 1:ncol(Contextual_data)){
  x <- cor.test(M_Census, Contextual_data[,i])
  Contextual_data_summary [i,] <- c(x$estimate, x$statistic, x$p.value)
}
row.names(Contextual_data_summary) <- colnames(Contextual_data)
Contextual_data_summary$AdjustedP_value <- p.adjust(Contextual_data_summary$p_value, method="fdr")

write.table(file="Tables/TableSContextualData_Summary.txt", x=Contextual_data_summary, sep="\t", quote=FALSE)

###########################################################
#### Thermophile Representation in CheckM Tree and IMG #### 
###########################################################
Finished_Mesophiles <- read.table("Input_Files/IMG_Finished_Mesophiles.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)
Finished_Thermophiles <- read.table("Input_Files/Finished_Thermophile_Genomes.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)

Genera <- unique(c(Finished_Mesophiles$Genus, Finished_Thermophiles$Genus))
Genus_Count_IMG <- data.frame(Genus = Genera, Thermophiles = rep(NA, length(Genera)), Mesophiles=rep(NA, length(Genera)))

for (i in 1:length(Genera)){
  Genus_Count_IMG$Thermophiles[i] <-sum(grepl(pattern = Genera[i], x= Finished_Thermophiles$Genus))
  Genus_Count_IMG$Mesophiles[i] <- sum(grepl(pattern = Genera[i], x= Finished_Mesophiles$Genus))
  }

GC_Subset <- Genus_Count_IMG[Genus_Count_IMG$Thermophiles>0,]
View(Genus_Count_IMG)

Families <- unique(c(Finished_Mesophiles$Family, Finished_Thermophiles$Family))
Family_Counts_IMG <- data.frame(Family=Families, Thermophiles=rep(NA, length(Families)), Mesophiles=rep(NA, length(Families)))
for (i in 1:length(Families)){
  Family_Counts_IMG$Thermophiles[i] <- sum(grepl(pattern = Families[i], x=Finished_Thermophiles$Family))
  Family_Counts_IMG$Mesophiles[i] <- sum(grepl(pattern=Families[i], x=Finished_Mesophiles$Family))
}


Orders <- unique(c(Finished_Mesophiles$Order, Finished_Thermophiles$Order))
Order_Counts_IMG <- data.frame(Order=Orders, Thermophiles=rep(NA, length(Orders)), Mesophiles=rep(NA, length(Orders)))
for (i in 1:length(Orders)){
  Order_Counts_IMG$Thermophiles[i] <- sum(grepl(pattern = Orders[i], x=Finished_Thermophiles$Order))
  Order_Counts_IMG$Mesophiles[i] <- sum(grepl(pattern=Orders[i], x=Finished_Mesophiles$Order))
}

Classes <- unique(c(Finished_Mesophiles$Class, Finished_Thermophiles$Class))
Class_Counts_IMG <- data.frame(Class=Classes, Thermophiles=rep(NA, length(Classes)), Mesophiles=rep(NA, length(Classes)))
for (i in 1:length(Classes)){
  Class_Counts_IMG$Thermophiles[i] <- sum(grepl(pattern = Classes[i], x=Finished_Thermophiles$Class))
  Class_Counts_IMG$Mesophiles[i] <- sum(grepl(pattern=Classes[i], x=Finished_Mesophiles$Class))
}


Phyla <- unique(c(Finished_Mesophiles$Phylum, Finished_Thermophiles$Phylum))
Phyla_Counts_IMG <- data.frame(Phylum=Phyla, Thermophiles=rep(NA, length(Phyla)), Mesophiles=rep(NA, length(Phyla)))
for (i in 1:length(Phyla)){
  Phyla_Counts_IMG$Thermophiles[i] <- sum(grepl(pattern = Phyla[i], x=Finished_Thermophiles$Phylum))
  Phyla_Counts_IMG$Mesophiles[i] <- sum(grepl(pattern=Phyla[i], x=Finished_Mesophiles$Phylum))
}


#################################################
#### Working with data for iTOL and MiGA Now ####
#################################################
Tree_Tips <- read.table("Input_Files/TreeTips.txt", sep="\t", header = FALSE, row.names = NULL, stringsAsFactors = FALSE)

Tree_Tips_1 <- Tree_Tips

colnames(Tree_Tips) <- "Genome_ID"
All_Archaea <- read.table("Input_Files/All_Archaea_IMG_April19_2018.txt",sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names=NULL)

All_Bacteria <- read.table("Input_Files/All_Bacteria_Test.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)

colnames(All_Archaea) <- colnames(All_Bacteria)

All_IMG <- rbind(All_Archaea, All_Bacteria)

colnames(All_IMG)[7] <- "Genome_ID"

All_IMG_Thermophiles <- All_IMG[All_IMG$Temperature.Range=="Thermophile",]
All_IMG_Thermophiles <- rbind(All_IMG_Thermophiles, All_IMG[All_IMG$Temperature.Range=="Hyperthermophile",])
All_IMG_Thermophiles <- rbind(All_IMG_Thermophiles, All_IMG[All_IMG$Temperature.Range=="Thermotolerant",])
All_IMG_Thermophiles_Finished <- All_IMG_Thermophiles[All_IMG_Thermophiles$Sequencing.Status=="Finished",]

Tree_Tips$Genome_ID <- gsub(pattern="IMG_", replacement = "", Tree_Tips$Genome_ID)

Joined_Tree_Data <- plyr::join(x=Tree_Tips, y=All_IMG, by="Genome_ID", type="left")

Complete_Joined_Tree_Data <- Joined_Tree_Data[complete.cases(Joined_Tree_Data$Sequencing.Status),]

MiGA_Results <- read.table("Input_Files/MiGA_MAG_Assignments.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

MAG_Taxonomy <- Complete_Clean$MarkerLineage
Unique_MAG_Taxonomy <- unique(MAG_Taxonomy)

MAG_Taxonomy_Phylum <- MiGA_Results$MiGA_Taxonomy_Simplified
MAG_Taxonomy_Phylum <- gsub("d__", "", MAG_Taxonomy_Phylum)
MAG_Taxonomy_Phylum <- gsub("p__", "", MAG_Taxonomy_Phylum)
MAG_Taxonomy_Phylum <- gsub("c__", "", MAG_Taxonomy_Phylum)
MiGA_Results$MAG[order(MiGA_Results$MAG)]==Complete_Clean$BinID[order(Complete_Clean$BinID)]

#Let's get rid of all those Pesky NA's in the dataset
Joined_Tree_Data$Phylum[is.na(Joined_Tree_Data$Phylum)] <- "None"
Joined_Tree_Data$SimplifiedTaxonomy <- Joined_Tree_Data$Phylum


Joined_Tree_Data$SimplifiedTaxonomy[Joined_Tree_Data$Phylum=="Proteobacteria"] <- Joined_Tree_Data$Class[Joined_Tree_Data$Phylum=="Proteobacteria"]

for (i in 1:nrow(Complete_Clean)){
  Joined_Tree_Data[Joined_Tree_Data$Genome_ID==MiGA_Results$MAG[i],][18] <- MAG_Taxonomy_Phylum[i]
}

Simplified_Tree_Data <- Joined_Tree_Data[,c(1,14,18)]

Simplified_Tree_Data$SimplifiedTaxonomy[is.na(Simplified_Tree_Data$SimplifiedTaxonomy)] <- "None"
length(unique(Simplified_Tree_Data$SimplifiedTaxonomy))



#Creat a color for every unique taxonomic rank
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_vector <- unique(col_vector)
 

Collapsed_Tree_Phyla <- read.table("Input_Files/Collapsed_Tree_Phyla.txt", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep="\t")

Collapsed_Tree_Phyla <- as.character(Collapsed_Tree_Phyla$V1)

Collapsed_Tree_Phyla <- Collapsed_Tree_Phyla[order(Collapsed_Tree_Phyla)]

Distinct_Colors <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080")



Reduced_Phyla <- Collapsed_Tree_Phyla[Collapsed_Tree_Phyla%in%MAG_Taxonomy_Phylum]




#Unique_Simplified_Taxa <- unique(Simplified_Tree_Data$SimplifiedTaxonomy)
#Unique_Simplified_Taxa <- Unique_Simplified_Taxa[order(Unique_Simplified_Taxa)]

Simplified_Tree_Data$Taxonomy_Colors <- rep("NA", nrow(Simplified_Tree_Data))

for(i in 1:length(Reduced_Phyla)){
  Simplified_Tree_Data$Taxonomy_Colors[Simplified_Tree_Data$SimplifiedTaxonomy==Reduced_Phyla[i]] <- Distinct_Colors[i]
}
Simplified_Tree_Data$Taxonomy_Colors[Simplified_Tree_Data$SimplifiedTaxonomy=="None"] <- "#FFFFFF"

Simplified_Tree_Data$Taxonomy_Colors[!Simplified_Tree_Data$SimplifiedTaxonomy%in%MAG_Taxonomy_Phylum] <- "#FFFFFF"

Tree_Taxonomy_Colors <- data.frame(Genome_ID=Tree_Tips_1$V1, TaxonomyColor=Simplified_Tree_Data$Taxonomy_Colors)


#Unique_Simplified_Taxa_Colors <- data.frame(Phylum=Collapsed_Tree_Phyla, TaxonomyColor=Taxonomy_Colors)



write.table(file="Figures/TreeInfo/Tree_Taxonomy_Colors.txt",x=Tree_Taxonomy_Colors, sep=",", col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.csv(file="Tree_Taxonomy.txt", x=matrix(Unique_Simplified_Taxa, nrow=1), quote=FALSE, row.names = FALSE)

#write.csv(file="Squares.txt", x=matrix(rep(1,length(Unique_Simplified_Taxa)), nrow=1), quote=FALSE, row.names=FALSE)

#write.csv(file="Legend_Colors.txt", x=matrix(Taxonomy_Colors, nrow=1),row.names=FALSE, quote=FALSE)

### Making Temperature Range Identifier for iTOL
unique(Joined_Tree_Data$Temperature.Range)
Simplified_Temp_Data <- Joined_Tree_Data
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="35 C", "Mesophile", Joined_Tree_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="30 C", "Mesophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="Mesophile 28C", "Mesophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="Psychrotrophic", "Psychrophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="Psychrotolerant", "Psychrophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="Hyperthermophile", "Thermophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="to 93 C", "Thermophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range <- gsub(pattern="Thermotolerant", "Thermophile", Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range[is.na(Simplified_Temp_Data$Temperature.Range)] <- "No Data"
unique(Simplified_Temp_Data$Temperature.Range)
Simplified_Temp_Data$Temperature.Range[Simplified_Temp_Data$Temperature.Range==unique(Simplified_Temp_Data$Temperature.Range)[3]] <- "No Data"
unique(Simplified_Temp_Data$Temperature.Range)

Temperature_Colors <- c("#FF0000", "#00FF00", "#FFFFFF", "#0000FF")
Temperature_Ranges <- unique(Simplified_Temp_Data$Temperature.Range)

Temperature_Color_Legend <- data.frame(Range=Temperature_Ranges, Color=Temperature_Colors)





Simplified_Temp_Tree_Data <- data.frame(Genome_ID=Tree_Tips_1$V1, Temperature_Range=Simplified_Temp_Data$Temperature.Range, Temp_Color= rep(NA, length(Tree_Tips_1$V1)))

for(i in 1:length(Temperature_Ranges)){
  Simplified_Temp_Tree_Data$Temp_Color[Simplified_Temp_Tree_Data$Temperature_Range==Temperature_Ranges[i]] <- Temperature_Colors[i] 
}

#Temp_Range_Legend_plot <- ggplot(Simplified_Temp_Tree_Data, aes(x=Temperature_Range))+
 # geom_bar(aes(fill=Temperature_Range), color="black")+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))+
  #scale_fill_manual(values=c("#00FF00", "#FFFFFF", "#0000FF","#FF0000"))
#ggsave(plot = Temp_Range_Legend_plot, filename = "Figures/Figures3A_Temp_Legend.eps", width = 6, height = 6, units = "in", device = "eps")

write.table(x=Simplified_Temp_Tree_Data[,c(1,3)], file = "Figures/TreeInfo/Tree_Temp_Range_Labels.txt", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE )


#### Make Color strips for MAG Identification on iTOL
MAG_Color <- rep("#8b0000", sum(grepl("METABAT", Tree_Tips_1$V1)))
MAG_Color[map_MG$Classification[apply(MAGS_Coverage[,1:12], 1, which.max)]!="FireAffected"] <- "#00468b"

MAG_Color <- data.frame(Tree_Tip=row.names(MAGS_Coverage), Tree_Tip_Color=MAG_Color)

IMG_Color <- data.frame(Tree_Tip=Tree_Tips_1$V1[!grepl(x=Tree_Tips_1$V1, pattern = "METABAT")], Tree_Tip_Color=rep("#ffffff", length(Tree_Tips_1$V1[!grepl(x=Tree_Tips_1$V1, pattern = "METABAT")])))

MAG_Color_Key <- data.frame(Tree_Tip_Color=c("#8b0000", "#00468b", "#ffffff"), Genome_Type=c("Fire Affected MAG", "Recovered Reference MAG", "IMG Genome"))

MAG_Identification<- rbind(IMG_Color, MAG_Color)

MAG_Color_Legend <- plyr::join(MAG_Color_Key, MAG_Identification, by="Tree_Tip_Color")

#MAG_Color_Legend_Plot <- ggplot(MAG_Color_Legend, aes(x=Genome_Type))+
 # geom_bar(aes(fill=Genome_Type), color="black")+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))+
  #scale_fill_manual(values=c("#8b0000", "#ffffff", "#00468b"))
#ggsave(plot=MAG_Color_Legend_Plot, file="Figures/Figures3A_GenomeType_Legend.eps", width=6, height = 6, units="in", device = "eps")




write.table(x=MAG_Identification, "Figures/TreeInfo/MAG_Identifying_Color.txt", quote = FALSE, col.names = FALSE, row.names = FALSE,sep=",")



MAGS_Coverage <- as.data.frame(MAGS_Coverage)
colnames(MAGS_Coverage) <- map_MG$Sample


PF_All_Centralia <- PF_All[PF_All$Phylum%in%Fixed_Phylum,]
#### Creating a dataset so we can show the distribution of bacteria, Archaea, and Proteobacteria
Additional_Bacteria <- PF_All_Centralia[PF_All_Centralia$Domain=="Bacteria",]
Additional_Bacteria$Phylum <- "Bacteria"

Additional_Archaea <- PF_All_Centralia[PF_All_Centralia$Domain=="Archaea",]
Additional_Archaea$Phylum <- "Archaea"

Genome_Size_IMG <- rbind(PF_All_Centralia, Additional_Bacteria, Additional_Archaea)




n <- as.data.frame(Phyla_Comm_Rel)
colnames(n) <- map_MG$Sample

n$Phylum <- row.names(n)


n.long <- melt(n, id.vars = "Phylum", measure.vars = colnames(n)[1:12])

n.long$Classification <- rep("FireAffected", nrow(n.long)) 
n.long$Classification[n.long$variable%in%map_MG$Sample[map_MG$Classification!="FireAffected"]] <- "Recovered/Reference"


n.long$Phylum <- factor(n.long$Phylum, levels= n$Phylum[order(Necessary_JGI[,9])])
n.long$Phylum <- factor(n.long$Phylum, levels = c(levels(n.long$Phylum)[c(31,33:35,1:30,32,36:41)],"Lentisphaerae"))

S2A <- ggplot(n.long, aes(x=Phylum, y=value, fill=Classification))+
  geom_boxplot()+
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.title.y = element_text(size=8))+
  labs(y="Relative Abundance")+
  scale_x_discrete(drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title.x = element_blank(),axis.title.y = element_text(size=8))+
  guides(fill=FALSE)
S2A

PF_All_Centralia <- PF_All[PF_All$Phylum%in%Fixed_Phylum,]
#### Creating a dataset so we can show the distribution of bacteria, Archaea, and Proteobacteria
Additional_Bacteria <- PF_All_Centralia[PF_All_Centralia$Domain=="Bacteria",]
Additional_Bacteria$Phylum <- "Bacteria"

Additional_Archaea <- PF_All_Centralia[PF_All_Centralia$Domain=="Archaea",]
Additional_Archaea$Phylum <- "Archaea"

Genome_Size_IMG_2 <- rbind(PF_All, Additional_Bacteria, Additional_Archaea)

Genome_Size_IMG_2 <- Genome_Size_IMG_2[Genome_Size_IMG_2$Phylum%in%levels(n.long$Phylum),]


Genome_Size_IMG_2$Phylum <- factor(Genome_Size_IMG_2$Phylum, levels = levels(n.long$Phylum))

S2B <- ggplot(Genome_Size_IMG_2, aes(x=Phylum, y=Genome.Size.....assembled/1000000))+
  geom_boxplot()+
  labs(y="Genome Size (Mbp)")+
  scale_x_discrete(drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title.x = element_blank(),axis.title.y = element_text(size=8))
S2B
library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(S2A), ggplotGrob(S2B), size="last"))

setEPS()
svg("Figures/FigureS2AB_16SPhylum.svg", width = 8, height=8)
par(ps = 12, cex = 1, cex.main = 1)
grid.draw(rbind(ggplotGrob(S2A), ggplotGrob(S2B), size="last"))
dev.off()

unique(MiGA_Results$MiGA.Taxonomy)[order(unique(MiGA_Results$MiGA.Taxonomy))]
Mag_Taxonomy_barplot <- MiGA_Results
Mag_Taxonomy_barplot$MiGA_Taxonomy_Simplified <- factor(Mag_Taxonomy_barplot$MiGA_Taxonomy_Simplified, levels=names(table(MiGA_Results$MiGA_Taxonomy_Simplified))[c(3:9, 11:15,1,2,10,16,17)])


Mag_Taxonomy_barplot$MAG_Type <- rep("FireAffected", nrow(Mag_Taxonomy_barplot))
Mag_Taxonomy_barplot$MAG_Type[map_MG$Classification[apply(MAGS_Coverage[,1:12], 1, which.max)]!="FireAffected"] <- "Recovered/Reference"


### Figure S2C 
#Mag_Taxonomy_barplot_plot <- ggplot(Mag_Taxonomy_barplot, aes(x=MiGA_Taxonomy_Simplified))+
 # geom_bar(aes(fill=MAG_Type))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title.x = element_blank(),axis.title.y = element_text(size=8))+
  #labs(y="Number of MAGs")+
  #geom_text(stat = "count",aes(label=..count.., x=MiGA_Taxonomy_Simplified),inherit.aes = FALSE,vjust=-.5, check_overlap = TRUE)+
  #scale_fill_manual(values = c("#8b0000", "#00468b"))
#ggsave(plot=Mag_Taxonomy_barplot_plot, file="Figures/FigureS2C_MAG_Taxonomy.svg", width=8, height = 5, units="in", device = "svg")
#ggsave(plot=Mag_Taxonomy_barplot_plot, file="Figures/FigureS2C_MAG_Taxonomy.eps", width=8, height = 5, units="in", device = "eps")


### RefSoil GenomeSize Boxplot
RefSoil_Genome_Size <- read.table("Input_Files/RefSoilGenomeSize.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
MAG_MiGA_Joined <- plyr::join(MAG_GenomeSize_Estimate, MiGA_Results, by="MAG")
MAG_MiGA_Joined$MiGA_Taxonomy_Simplified <- gsub("d__", "", MAG_MiGA_Joined$MiGA_Taxonomy_Simplified)
MAG_MiGA_Joined$MiGA_Taxonomy_Simplified <- gsub("p__", "", MAG_MiGA_Joined$MiGA_Taxonomy_Simplified)
MAG_MiGA_Joined$MiGA_Taxonomy_Simplified <- gsub("c__", "", MAG_MiGA_Joined$MiGA_Taxonomy_Simplified)

RefSoil_Genome_Size$SimplifiedPhylum <- RefSoil_Genome_Size$Phylum
RefSoil_Genome_Size$SimplifiedPhylum[RefSoil_Genome_Size$SimplifiedPhylum=="Proteobacteria"] <- RefSoil_Genome_Size$Class[RefSoil_Genome_Size$SimplifiedPhylum=="Proteobacteria"]
unique(RefSoil_Genome_Size$SimplifiedPhylum)

temp_proteo <- RefSoil_Genome_Size[RefSoil_Genome_Size$Phylum=="Proteobacteria", ]
temp_proteo$SimplifiedPhylum <- temp_proteo$Phylum

temp_bac <- RefSoil_Genome_Size[!grepl("archaeota", RefSoil_Genome_Size$Phylum),]
temp_bac$SimplifiedPhylum<- "Bacteria"

temp_arc <- RefSoil_Genome_Size[grepl("archaeota", RefSoil_Genome_Size$Phylum),]
temp_arc$SimplifiedPhylum <- "Archaea"

RefSoil_Genome_Size <- rbind(RefSoil_Genome_Size, temp_proteo, temp_bac, temp_arc)

Reduced_RefSoil_Genome_Size <- RefSoil_Genome_Size[RefSoil_Genome_Size$SimplifiedPhylum%in%MAG_MiGA_Joined$MiGA_Taxonomy_Simplified,]
Reduced_RefSoil_Genome_Size <- rbind(Reduced_RefSoil_Genome_Size, temp_arc)

table(Reduced_RefSoil_Genome_Size$SimplifiedPhylum)

"#8b0000", "#00468b"
MAG_Color <- rep("#8b0000", nrow(MAG_MiGA_Joined))
MAG_Color[map_MG$Classification[apply(MAGS_Coverage[,1:12], 1, which.max)]!="FireAffected"] <- "#ffffff"

MAG_MiGA_Joined <- MAG_MiGA_Joined[order(MAG_MiGA_Joined$MAG),]

MAG_MiGA_Joined$MAG==row.names(MAGS_Coverage)

 

Reduced_RefSoil_Genome_Size$SimplifiedPhylum <- factor(Reduced_RefSoil_Genome_Size$SimplifiedPhylum, levels = c("Bacteria", "Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Chlamydiae", "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Alphaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Archaea","Crenarchaeota","Thaumarchaeota"))

Reduced_RefSoil_Genome_Size$Redundant <- rep("No", nrow(Reduced_RefSoil_Genome_Size))
Reduced_RefSoil_Genome_Size$Redundant[Reduced_RefSoil_Genome_Size$SimplifiedPhylum%in%c("Proteobacteria", "Bacteria", "Archaea")] <- "Yes"


MAG_MiGA_Joined$MAGClass <- rep("FireAffected", nrow(MAG_MiGA_Joined))
MAG_MiGA_Joined$MAGClass[map_MG$Classification[apply(MAGS_Coverage[,1:12], 1, which.max)]!="FireAffected"] <- "Recovered/Reference"


MAG_MiGA_Joined$MiGA_Taxonomy_Simplified <- factor(MAG_MiGA_Joined$MiGA_Taxonomy_Simplified, levels = c("Bacteria", "Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Chlamydiae", "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Alphaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Archaea","Crenarchaeota","Thaumarchaeota"))



MAG_MiGA_Joined_Reduced <- MAG_MiGA_Joined[,c(1,2,6,7)]
Reduced_RefSoil_Genome_Size_Reduced <- Reduced_RefSoil_Genome_Size[,c(1,10,11)]
Reduced_RefSoil_Genome_Size_Reduced$MAGClass <- rep(NA, nrow(Reduced_RefSoil_Genome_Size_Reduced))

colnames(Reduced_RefSoil_Genome_Size_Reduced) <- colnames(MAG_MiGA_Joined_Reduced)

Figure3b_Data <- rbind(Reduced_RefSoil_Genome_Size_Reduced, MAG_MiGA_Joined_Reduced)

Figure3b_Data$MiGA_Taxonomy_Simplified <- factor(Figure3b_Data$MiGA_Taxonomy_Simplified, levels = c("Bacteria", "Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Chlamydiae", "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Alphaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Archaea","Crenarchaeota","Thaumarchaeota"))

Figure3b_Data$Redundant <- rep("No", nrow(Figure3b_Data))
Figure3b_Data$Redundant[grepl("RefSoil",Figure3b_Data$MAG)&Figure3b_Data$MiGA_Taxonomy_Simplified%in%c("Bacteria", "Archaea", "Proteobacteria")] <- "Yes"



Figure3B_Plot <- ggplot(Figure3b_Data, aes(x=MiGA_Taxonomy_Simplified, y=Genome_Length/1000000))+
  geom_boxplot(data=Figure3b_Data[grep("RefSoil", Figure3b_Data$MAG),], aes(fill=Redundant))+
  geom_jitter(data=Figure3b_Data[grep("METABAT", Figure3b_Data$MAG),], shape=2, aes(x=MiGA_Taxonomy_Simplified, y=Genome_Length/1000000, color=MAGClass))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title.x = element_blank(),axis.title.y = element_text(size=8))+
  scale_fill_manual(values=c("#FFFFFF","#A9A9A9"))+
  scale_color_manual(values=c("#8b0000", "#00468b"))+
  scale_y_continuous(limits = c(0,15))+
  labs(y="Genome Size (Mbp)")+
  scale_x_discrete(breaks=c("Bacteria", "Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Chlamydiae", "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Alphaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Archaea","Crenarchaeota","Thaumarchaeota"), drop=FALSE)+
  geom_text(data=Figure3b_Data[grep("RefSoil", Figure3b_Data$MAG),], stat = "count",aes(label=..count.., x=MiGA_Taxonomy_Simplified),inherit.aes = FALSE,vjust=-.5, check_overlap = TRUE, y=1)
ggsave(plot = Figure3B_Plot, file="Figures/Figure3B_OLD_GenomeBoxplot.eps", )

MAG_Bac_Dupes <- Figure3b_Data[grepl("METABAT", Figure3b_Data$MAG)&Figure3b_Data$MiGA_Taxonomy_Simplified%in%c("Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Chlamydiae", "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Alphaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia"), ]

MAG_Bac_Dupes$MiGA_Taxonomy_Simplified <- "Bacteria"

MAG_Arc_Dupes <- Figure3b_Data[grepl("METABAT", Figure3b_Data$MAG)&Figure3b_Data$MiGA_Taxonomy_Simplified%in%c("Crenarchaeota","Thaumarchaeota"), ]
MAG_Arc_Dupes$MiGA_Taxonomy_Simplified <- "Archaea"


MAG_Proteo_Dupes <- Figure3b_Data[grepl("METABAT", Figure3b_Data$MAG)&Figure3b_Data$MiGA_Taxonomy_Simplified%in%c("Alphaproteobacteria", "Gammaproteobacteria"), ]
MAG_Proteo_Dupes$MiGA_Taxonomy_Simplified <- "Proteobacteria"


Figure3b_New_Data <- rbind(Figure3b_Data, MAG_Bac_Dupes, MAG_Arc_Dupes, MAG_Proteo_Dupes)

Figure3b_New_Data$MAGClass[is.na(Figure3b_New_Data$MAGClass)] <- "RefSoil"

Figure3b_New_Data$MAGClass[Figure3b_New_Data$MAGClass=="Recovered/Reference"] <- "Ambient"

Figure3b_New_Data$MAGClass <- factor(Figure3b_New_Data$MAGClass, levels=c("RefSoil", "FireAffected", "Ambient"))


BoxplotColors <- data.frame(Phylum=unique(MAG_Taxonomy_Phylum)[order(unique(MAG_Taxonomy_Phylum))], Color=Distinct_Colors[1:17])
#BoxplotColors <- Taxonomy_Legend[Taxonomy_Legend$Phylum%in%Figure3b_New_Data$MiGA_Taxonomy_Simplified,]
BoxplotColors$Phylum <- as.character(BoxplotColors$Phylum)
BoxplotColors$Color <- as.character(BoxplotColors$Color)

BoxplotColors <- rbind(BoxplotColors, c("Archaea", "#808080"))

BoxplotColors_2 <- BoxplotColors$Color
names(BoxplotColors_2) <- BoxplotColors$Phylum

Figure3b_New_plot
Figure3b_New_Freey_plot <- ggplot(Figure3b_New_Data, aes(x=MAGClass, y=Genome_Length/1000000, fill=MiGA_Taxonomy_Simplified))+
  geom_boxplot(outlier.size=.5, outlier.shape = 1)+
  geom_jitter(data=Figure3b_New_Data[grepl("METABAT",Figure3b_New_Data$MAG),], aes(x=MAGClass, y=Genome_Length/1000000), width = .25, size=.5 )+
  facet_wrap(~MiGA_Taxonomy_Simplified, nrow = 6, scales="free_y")+
  theme_bw(base_size = 7)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title.x = element_blank(),axis.title.y = element_text(size=8), panel.grid.minor = element_blank())+
  labs(y="Genome Size (Mbp)")+
  scale_fill_manual(values=BoxplotColors_2)#+
  #scale_y_continuous(limits = c(0,15))+
  #geom_text(data=Figure3b_New_Data[grep("RefSoil", Figure3b_New_Data$MAGClass),], stat = "count",aes(label=..count.., x=MAGClass),vjust=-.5, check_overlap = TRUE, y=2, size=3)+
  #geom_text(data=Figure3b_New_Data[grep("FireAffected", Figure3b_New_Data$MAGClass),], stat = "count",aes(label=..count.., x=MAGClass),vjust=-.5, check_overlap = TRUE, y=2, size=3)+
  #geom_text(data=Figure3b_New_Data[grep("Ambient", Figure3b_New_Data$MAGClass),], stat = "count",aes(label=..count.., x=MAGClass),vjust=-.5, check_overlap = TRUE, y=2, size=3)
Figure3b_New_Freey_plot
ggsave("Figures/Figure3b_Free_y_Faceted_Boxplots.eps", Figure3b_New_Freey_plot, height = 5, width=5, units = "in", device="eps")

#Figure3b_New_Fixed_y_plot <- ggplot(Figure3b_New_Data, aes(x=MAGClass, y=Genome_Length/1000000, fill=MiGA_Taxonomy_Simplified))+
 # geom_boxplot()+
  #geom_jitter(data=Figure3b_New_Data[grepl("METABAT",Figure3b_New_Data$MAG),], aes(x=MAGClass, y=Genome_Length/1000000) )+
  #facet_wrap(~MiGA_Taxonomy_Simplified, nrow = 6)+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title.x = element_blank(),axis.title.y = element_text(size=8))+
  #labs(y="Genome Size (Mbp)")+
  #scale_fill_manual(values=BoxplotColors_2)+
  #scale_y_continuous(limits = c(0,15))+
  #geom_text(data=Figure3b_New_Data[grep("RefSoil", Figure3b_New_Data$MAGClass),], stat = "count",aes(label=..count.., x=MAGClass),vjust=-.5, check_overlap = TRUE, y=2, size=6)+
  #geom_text(data=Figure3b_New_Data[grep("FireAffected", Figure3b_New_Data$MAGClass),], stat = "count",aes(label=..count.., x=MAGClass),vjust=-.5, check_overlap = TRUE, y=2, size=6)+
  #geom_text(data=Figure3b_New_Data[grep("Ambient", Figure3b_New_Data$MAGClass),], stat = "count",aes(label=..count.., x=MAGClass),vjust=-.5, check_overlap = TRUE, y=2, size=6)
#ggsave("Figures/Figure3b_Fixed_y_Faceted_Boxplots.eps", Figure3b_New_Fixed_y_plot, height = 11, units = "in", device="eps")







###
colnames(MiGA_Results)
colnames(Complete_Clean)
colnames(Complete_Clean)[1] <- "MAG"
TableS6 <- plyr::join(Complete_Clean, MiGA_Results, by="MAG")

TableS6 <- TableS6[,c(1,12,13,15)]
write.table(x=TableS6, file = "Tables/TableS6_MAG_Assignments.txt", sep="\t", col.names = TRUE, row.names = FALSE,quote = FALSE)


####
#### MAGs coverage for iTOL ####

C10_Max_MAGS <- MAGS_Coverage[apply(MAGS_Coverage,1,which.max)==7,]
row.names(C10_Max_MAGS)


MAGS_Coverage_iTOL <- decostand(MAGS_Coverage, method = "total", MARGIN = 1) 
rowSums(MAGS_Coverage_iTOL)

Tree_Tips_TempDesignations <- Tree_Tips_1[grep("IMG_", Tree_Tips_1$V1), 1]

IMG_Temp_Matrix <- matrix(data=0,nrow=length(Tree_Tips_TempDesignations), ncol=16)
IMG_Temp_Matrix <- as.data.frame(IMG_Temp_Matrix)
colnames(IMG_Temp_Matrix) <- c(colnames(MAGS_Coverage_iTOL), "Thermophile", "Mesophile", "Psychrophile","No Data")
row.names(IMG_Temp_Matrix) <- Tree_Tips_TempDesignations                               

x <- Simplified_Temp_Data[!grepl("META", Simplified_Temp_Data$Genome_ID),]

IMG_Temp_Matrix$Thermophile[x$Temperature.Range=="Thermophile"] <- 1
IMG_Temp_Matrix$Mesophile[x$Temperature.Range=="Mesophile"] <- 1
IMG_Temp_Matrix$Psychrophile[x$Temperature.Range=="Psychrophile"] <-1  
IMG_Temp_Matrix$`No Data`[x$Temperature.Range=="No Data"] <- 1

MAGS_Coverage_iTOL$Thermophile <- 0
MAGS_Coverage_iTOL$Mesophile <- 0  
MAGS_Coverage_iTOL$Psychrophile <- 0
MAGS_Coverage_iTOL$`No Data` <- 0

Tree_Abundance_Matrix <- rbind(MAGS_Coverage_iTOL, IMG_Temp_Matrix)

write.table(x=MAGS_Coverage_iTOL, file="Figures/TreeInfo/Tree_MAG_Abundance.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(x=IMG_Temp_Matrix, file="Figures/TreeInfo/Tree_IMG_Temp.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(x=Tree_Abundance_Matrix, file = "Figures/TreeInfo/Tree_Abundance_Matrix.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep=",")


paste("<h1>","Homo sapiens","</h1>", sep="")

Simplified_Tree_Data$Genome_IDpaste("<h1>",Simplified_Tree_Data$SimplifiedTaxonomy,"</h1>", sep="")

Hover_dataframe <- data.frame(Node=Tree_Tips_1$V1, Phylum="Phylum",Info=paste("<h1>",Simplified_Tree_Data$SimplifiedTaxonomy,"</h1>","<h2>",Tree_Tips_1$V1,"</h2>", sep=""))
write.table(x=Hover_dataframe, file = "Figures/TreeInfo/Hover_Info.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep=",")

##############################
#### Database Comparisons ####
##############################
# JGI Enzyme, pfam, COG datatables
enzyme <- read.table("Input_Files/abundance_enzyme_79244.txt", sep="\t", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
colnames(enzyme) <- map_MG$Sample
enzyme.rel <- decostand(enzyme, MARGIN = 2, method="total")
enzyme.d <- vegdist(t(enzyme.rel), method = "bray")

cog <- read.table("Input_Files/abundance_cog_79133.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
colnames(cog) <- map_MG$Sample
cog.rel <- decostand(cog, MARGIN=2, method="total")
cog.d<- vegdist(t(cog.rel), method="bray")

pfam <- read.table("Input_Files/abundance_pfam_72927.txt", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
colnames(pfam) <- map_MG$Sample
pfam.rel <- decostand(pfam, MARGIN=2, method="total")
pfam.d <- vegdist(t(pfam.rel), method="bray")

Mod.rel <- decostand(mid_mod, MARGIN=2, method="total")
Mod.d <- vegdist(t(Mod.rel), method="bray")

KO.d <- vegdist(t(KO.rel), method="bray")

#Correlation between Functional Potential matrices 
M_v_cog.mantel <- mantel(Mod.d, cog.d)
M_v_enzyme.mantel <- mantel(Mod.d, enzyme.d)
M_v_pfam.mantel <- mantel(Mod.d, pfam.d)
cog_v_enzyme.mantel <- mantel(cog.d, enzyme.d)
cog_v_pfam.mantel <- mantel(cog.d, pfam.d)
enzyme_v_pfam.mantel <- mantel(enzyme.d, pfam.d)
KO_v_cog.mantel <- mantel(KO.d, cog.d)
KO_v_enzyme.mantel <- mantel(KO.d, enzyme.d)
KO_v_pfam.mantel <- mantel(KO.d, pfam.d)

#Distance1 <-c("Module(BrayCurtis)", "Module(BrayCurtis)", "Module(BrayCurtis)", "KO(BrayCurtis)", "KO(BrayCurtis)", "UniFrac")
#Distance2<- c("KO(BrayCurtis)", "UniFrac", "Space","UniFrac", "Space", "Space")
#Mantel_Names<- c("Distance1", "Distance2","Mantel_R", "p-value")
#Mantel_R <- c(M_v_K.mantel[3], M_v_U.mantel[3], M_v_Space.mantel[3], K_v_U.mantel[3], K_v_Space.mantel[3], U_v_Space.mantel[3])
#Mantel_Pvalues <- c(M_v_K.mantel[4], M_v_U.mantel[4], M_v_Space.mantel[4], K_v_U.mantel[4], K_v_Space.mantel[4], U_v_Space.mantel[4])
#Mantel_Summary <- cbind(Distance1, Distance2, Mantel_R, Mantel_Pvalues)

pfam <- read.table("Input_Files/pfam_abund_24July2018.txt", sep="\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
pfam


pfam_ribosomal <- read.table("Input_Files/pfam_ribosomal_families.txt", sep="\t", header=FALSE, row.names=NULL, stringsAsFactors = FALSE)

pfam.rel <- decostand(x=pfam, method="total", MARGIN=2)

pfam_scg_Absolute <- pfam[pfam_ribosomal$V1,]
pfam_scg_rel <- pfam.rel[pfam_ribosomal$V1,]

# For use in Odds Ratios
Average_MG_SCG.pfam <- apply(pfam_scg_rel, 1, mean)
#For use in Normalizing KEGG Data
Average_SCG.pfam <- apply(pfam_scg_Absolute, 2, median)

pfam.sr<- NULL
for(i in 1:nrow(pfam)){
  pfam.sr <- rbind(pfam.sr, pfam[i,]/Average_SCG.pfam)
}

plasmid_pfams <- read.table("~/Downloads/pfam_Plamsids_Jorgensenetal_plosone.txt", sep="\t", header=FALSE, row.names=NULL, stringsAsFactors = FALSE)

pfam.sr.plasmid <- pfam.sr[plasmid_pfams$V1,]
pfam.sr.plasmid <- pfam.sr.plasmid[rowSums(pfam.sr.plasmid)>0,]
colnames(pfam.sr.plasmid) <- map_MG$Sample

zed <- cor.test(y=as.numeric(pfam.sr.plasmid[1,]), x=map_MG$SoilTemperature_to10cm)
zed$statistic
zed$estimate

plasmid_cor_test_Results <- data.frame(Pfam=row.names(pfam.sr.plasmid), Estimate=rep(NA, 11), Pvalue= rep(NA,11))

for (i in 1:nrow(pfam.sr.plasmid)){
  zed <- cor.test(y=as.numeric(pfam.sr.plasmid[i,]), x=map_MG$SoilTemperature_to10cm)
  plasmid_cor_test_Results$Estimate[i] <- zed$estimate
  plasmid_cor_test_Results$Pvalue[i] <- zed$p.value
}

write.table(x=plasmid_cor_test_Results, file="Tables/Plasmid_PFAM_Correlations.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

KO <- as.data.frame(KO)
Odds_Ratio_pfam <- NULL
for (i in 1:nrow(pfam_scg_rel)){
  Odds_Ratio_pfam <- rbind(Odds_Ratio_pfam, pfam_scg_rel[row.names(pfam_scg_rel[i,]),]/Average_MG_SCG.pfam[i])
}

scg_pfam_OR_Correlations <- data.frame(Pfam=row.names(Odds_Ratio_pfam), PearsonsR=rep(NA, nrow(Odds_Ratio_pfam)), Pvalue=rep(NA, nrow(Odds_Ratio_pfam)))
for(i in 1:nrow(Odds_Ratio_pfam)){
  bravo <- cor.test(y=as.numeric(Odds_Ratio_pfam[i,]), x=map_MG$SoilTemperature_to10cm)
  scg_pfam_OR_Correlations$PearsonsR[i] <- bravo$estimate
  scg_pfam_OR_Correlations$Pvalue[i]<- bravo$p.value
}


######################
### Boxplot Tables ###
######################
PF_All_2 <- rbind(PF_All, Additional_Archaea, Additional_Bacteria)

Phylum_For_BP_Table <- c(row.names(Phyla_Comm_Rel), "Lentisphaerae", "cadidate division SR1")

BP_Table <- data.frame(Phylum=Phylum_For_BP_Table, Genomes=rep(0, length(Phylum_For_BP_Table)))
for(i in 1:length(Phylum_For_BP_Table)){
  BP_Table$Genomes[i] <- length(PF_All_2$Phylum[PF_All_2$Phylum==Phylum_For_BP_Table[i]])
}

write.table(file = "Tables/TableSX_PhylumGenomeNumbers.txt", x=BP_Table, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

