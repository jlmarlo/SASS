# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:38:06 2023

@author: jillian
"""
import argparse


parser = argparse.ArgumentParser(description='Script to create slurm scripts for creation of simulated genomes')
parser.add_argument('settings',type=str,help='tab separated file including settings for the generation of genomes')
parser.add_argument('reference',type=str,help='Location of reference file')
parser.add_argument('email',type=str,help='Email address where slurm progress alerts will be sent')
parser.add_argument('shortPart',type=str,help='job partition that will be used for genomes that take < 4 days')
parser.add_argument('longPart',type=str,help='job partition that will be used for genomes that take > 4 days')

args = parser.parse_args()

#Read in generation settings
settingsData = open(args.settings,'r')
logFileDir = ".logs/Pipeline/"
referencelocation = args.reference
email = args.email
shortPart = args.shortPart## to be used for jobs taking less than 4 days
longPart = args.longPart ## to be used for jobs taking longer than 4 days

#Loop through each line
header = settingsData.readline()
for line in settingsData:
    line = line.rstrip()
#Remove enter
    emp,seed1,coverage,frag = line.split('\t')
#Split data into categories
#Write slurm script 
    for x in range(2):
        if x ==0:
            seed = seed1
        else:    
            seed = seed2
        output1name = "Seed" +seed+"Generating.slurm"
        script = open(output1name,"w") 
        if int(coverage) > 17:
            script.write(
                "#!/bin/bash -l \n"
                + "#SBATCH -t 168:00:00 \n"
                + "#SBATCH -p "+longPart+" \n")
        else:
            script.write(
                "#!/bin/bash -l \n"
                + "#SBATCH -t 96:00:00 \n"
                + "#SBATCH -p "+shortPart+" \n")
        script.write("#SBATCH --nodes=1 \n"
                + "#SBATCH --ntasks=1 \n"
                + "#SBATCH --cpus-per-task=1 \n"
                + "#SBATCH --mem=60gb \n"
                + "#SBATCH --mail-type=ALL \n"
                + "#SBATCH --mail-user="+email+"\n"
                + "#SBATCH --job-name GeneratingSeed" + seed + "\n"
                + "#SBATCH -o "+logFileDir+"%j.Seed"+seed + ".out\n"
                + "#SBATCH -e "+logFileDir+"%j.Seed" +seed + ".err\n"
                + "\n"
                + "conda activate Python2\n"
                + "python neat_genreads/utilities/genReads.py \\\n"
                + "\t-r"+referencelocation +" \\\n"
                + "\t-R " +frag +" \\\n"
                + "\t-o outputs/Genomes/Seed" +seed+ " \\\n"
                + "\t-c "+coverage+ " \\\n"
                + "\t-e outputs/Models/ErrorModels/" +horse+"ErrorModel \\\n"
                + "\t-M 0 \\\n"
                + "\t-v outputs/SeedVCFs/Seed" +seed+"Golden.vcf \\\n"
                + "\t--pe-model outputs/Models/FragLengthModels/"+horse+"/fraglen.p \\\n"
                + "\t--gc-model outputs/Models/GCModels/"+horse+"GCModel \\\n"
                + "\t--vcf \\\n"
                +"\t--bam \\\n"
                + "\t --gz"
                + "\n\n")

