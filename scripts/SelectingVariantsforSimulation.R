args <- commandArgs(trailingOnly=TRUE)
##arg 1 - path to director containing necessary libraries
##arg 2 - path to file with snps
##arg 3 - path to file with indels
##arg 4 - number of snps to insert
##arg 5 - number of indels to insert
##arg 6 - output prefix
##arg 7 - number of samples to select for
##Arg 8 - number from 1- infinite to set initial seed

#args<- c('~/R-library','VariantSelectionCodes/VariantsfromRob9M.tsv','VariantSelectionCodes/IndelsforSimulation.tsv',4000000,400000,'TestRun',10,243581)
.libPaths(args[1])
# List of required packages
packages <- c("ggplot2", "tidyverse")

# Function to check if a package is installed, and if not, install it
check_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Apply the function to each package in the list
lapply(packages, check_install)

# Confirm installation
cat("All packages are installed and loaded.\n")

variants <- read_tsv(args[2])
indels <- read_tsv(args[3])

#Combine all unknowns
variants$CHROMCn <- variants$CHROM
variants$CHROMCn <- gsub(".*chrUn.*", "chrUn",variants$CHROM)
#Remove Unknowns because the pipeline won't find them anyway
variants <- variants %>% filter(CHROMCn!='chrUn')
totalvariants <- nrow(variants)

table(variants$CHROMCn)

##Combine unknown indels

indels$CHROMCn <- indels$CHROM
indels$CHROMCn <- gsub(".*chrUn.*", "chrUn",indels$CHROM)
indels <- indels %>% filter(CHROMCn!='chrUn')
indels <- indels %>% filter(CHROMCn!='MSY')
indels <- indels %>% filter(CHROMCn!='chrM')
totalindels <- nrow(indels)


#Creating proportions and expected numbers
numbers<-data.frame(table(variants$CHROMCn))
numbers$proportions <- (numbers$Freq/totalvariants)
numbers$target <- numbers$proportions*(as.integer(args[4]))
numbers$target <- round(numbers$target,digits=0)


snpName <- paste("./",args[6],'SNPSpecs.tsv',sep='')

write.table(numbers,snpName,quote=FALSE,row.names=FALSE)

#Creating proportions and expected numbers for indels
numbersindels<-data.frame(table(indels$CHROMCn))
numbersindels$proportions <- (numbersindels$Freq/totalindels)
numbersindels$target <- numbersindels$proportions*(as.integer(args[5]))
numbersindels$target <- round(numbersindels$target,digits=0)

IndeName <- paste("./",args[6],'IndelSpecs.tsv',sep='')

write.table(numbersindels,IndeName,quote=FALSE,row.names=FALSE)





# Set the seed
set.seed(args[8])
print(paste("This is the seed used to generate seeds for randomized selection",as.integer(args[8]))) ##setting seed to randomly select seeds
##change this so that a seed is randomly generated but then it is printed

#creating list of random seeds
seeds <- sample(100:999, args[7])

##update to include prefix
write.table(seeds,"./SeedList.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#randomly selection snps and indels to meet criteria for each seed
samplevariants<- data.frame() #create snp dataframe
samplevariantsindels <- data.frame() #create indel dataframe
samplevariantsboth <-data.frame() #create dataframe for both snps and indels
for (seed in seeds){ #for every seed/simulated genome
  set.seed(seed) #set the seed for variant selection
  for (chr in numbers$Var1){ 
    TargetVariants <- variants %>% filter(CHROMCn ==chr) #create list of variants only on desired chromosome
    selectedvariants<- TargetVariants[sample(nrow(TargetVariants),size=numbers$target[which(numbers$Var1==chr)]),] #randomly sample variants until meet target number of variants for that chromosome
    samplevariants <- rbind(samplevariants,selectedvariants) ##add to list of variants for the seed/simulated genome
  }
  
  filename <- paste("./VariantLists/seed", seed,'variantregionlistSNPs.tsv', sep='') #save
  write.table(samplevariants[,1:2],filename,quote=FALSE,row.names=FALSE,sep='\t',col.names=FALSE) #save
  set.seed(seed) #reset seed for indel selection
  for (chr in numbersindels$Var1){
    TargetVariants <- indels %>% filter(CHROMCn ==chr)
    selectedvariants<- TargetVariants[sample(nrow(TargetVariants),size=numbersindels$target[which(numbersindels$Var1==chr)]),]
    samplevariantsindels <- rbind(samplevariantsindels,selectedvariants)
  }
  filename <- paste("./VariantLists/seed", seed,'variantregionlistindels.tsv', sep='') #save indels
  write.table(samplevariantsindels[,1:2],filename,quote=FALSE,row.names=FALSE,sep='\t',col.names=FALSE) #save indels
  samplevariantsboth <- rbind(samplevariants,samplevariantsindels) # combine snps and indels
  filename <- paste("./Lists/seed", seed,'variantlist.tsv', sep='')#save snps and indels
  write.table(samplevariantsboth[,1:4],filename,quote=FALSE,row.names=FALSE,sep='\t') #save both
  samplevariants <- data.frame() #clear frames for next seed
  samplevariantsindels <- data.frame()
  samplevariantsboth<-data.frame()
}
