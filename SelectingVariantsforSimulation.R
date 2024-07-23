args <- commandArgs(trailingOnly=TRUE)

SNPfile <- args[1]
Indels <- args[2]

.libPaths('~/R-library')
library(tidyverse)
library(psych)
library(dplyr)


variants <- read_tsv(SNPfile)
indels <- read_tsv(Indels)

#Combine all unknowns
variants$CHROMCn <- variants$CHROM
variants$CHROMCn <- gsub(".*chrUn.*", "chrUn",variants$CHROM)

#Remove Unknowns because the pipeline won't find them anyway
variants <- variants %>% filter(CHROMCn!='chrUn')
totalvariants <- nrow(variants)

##Combine unknown indels

indels$CHROMCn <- indels$CHROM #create new chromosome column
indels$CHROMCn <- gsub(".*chrUn.*", "chrUn",indels$CHROM) ##combine unknown
indels <- indels %>% filter(CHROMCn!='chrUn') #remove unknown
indels <- indels %>% filter(CHROMCn!='MSY') #remove MSY (no MSY in SNPs)
indels <- indels %>% filter(CHROMCn!='chrM') #remove mitochondrial
totalindels <- nrow(indels)


#Creating proportions and expected numbers
numbers<-data.frame(table(variants$CHROMCn))
numbers$proportions <- (numbers$Freq/totalvariants)
numbers$target <- numbers$proportions*4000000
numbers$target <- round(numbers$target,digits=0)

write.table(numbers,"~/Minn 2021-2022/SimulatingData90August/NumberSpecificationsRob.tsv")

#Creating proportions and expected numbers for indels
numbersindels<-data.frame(table(indels$CHROMCn))
numbersindels$proportions <- (numbersindels$Freq/totalindels)
numbersindels$target <- numbersindels$proportions*400000
numbersindels$target <- round(numbersindels$target,digits=0)

write.table(numbersindels,"~/Minn 2021-2022/SimulatingData90August/IndelNumberSpecifications.tsv")

set.seed(4643)

#creating seeds for all 15 samples
#For third round of simulations. (Doing 90 genomes) changed from 1:500 to 500:999
seeds <- sample(1:1000000, 90)
write.table(seeds,"~/Minn 2021-2022/SimulatingData90August/Seedlist90Genomes.tsv",row.names=FALSE, col.names=FALSE)


#Choosing correct number of snps from the correct chromosome.
samplevariants<- data.frame()
samplevariantsindels <- data.frame()
samplevariantsboth <-data.frame()


for (seed in seeds){
  set.seed(seed)
  for (chr in numbers$Var1){
    TargetVariants <- variants %>% filter(CHROMCn ==chr)
    selectedvariants<- TargetVariants[sample(nrow(TargetVariants),size=numbers$target[which(numbers$Var1==chr)]),]
    samplevariants <- rbind(samplevariants,selectedvariants)
  }
  
  filename <- paste("~/Minn 2021-2022/SimulatingData90August/seed", seed,'variantregionlistSNPs.tsv', sep='')
  write.table(samplevariants[,1:2],filename,quote=FALSE,row.names=FALSE,sep='\t',col.names=FALSE)
  set.seed(seed)
  for (chr in numbersindels$Var1){
    TargetVariants <- indels %>% filter(CHROMCn ==chr)
    selectedvariants<- TargetVariants[sample(nrow(TargetVariants),size=numbersindels$target[which(numbersindels$Var1==chr)]),]
    samplevariantsindels <- rbind(samplevariantsindels,selectedvariants)
  }
  filename <- paste("~/Minn 2021-2022/SimulatingData90August/seed", seed,'variantregionlistindels.tsv', sep='')
  write.table(samplevariantsindels[,1:2],filename,quote=FALSE,row.names=FALSE,sep='\t',col.names=FALSE)
  samplevariantsboth <- rbind(samplevariants,samplevariantsindels)
  filename <- paste("~/Minn 2021-2022/SimulatingData90August/seed", seed,'variantlist.tsv', sep='')
  write.table(samplevariantsboth[,1:4],filename,quote=FALSE,row.names=FALSE,sep='\t')
  samplevariants <- data.frame()
  samplevariantsindels <- data.frame()
  samplevariantsboth<-data.frame()
}
