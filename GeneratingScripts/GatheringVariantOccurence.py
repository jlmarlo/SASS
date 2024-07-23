# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:38:06 2023

@author: jillian
"""
#Read in data
settingsData = open('GenerationSettings.txt','r')
#Loop through each line
header = settingsData.readline()
for line in settingsData:
    line = line.rstrip()
#Remove enter
    horse,seed1,seed2,coverage,frag = line.split('\t')
#Split data into categories
#Write slurm script 
    for x in range(2):
        if x ==0:
            seed = seed1
        else:    
            seed = seed2
        output1name = "Seed" +seed+"Generating.slurm"
        script = open(output1name,"w") 
        if coverage == '17':
            script.write(
                "#!/bin/bash -l \n"
                + "#SBATCH -t 168:00:00 \n"
                + "#SBATCH -p max\n")
        else:
            script.write(
                "#!/bin/bash -l \n"
                + "#SBATCH -t 96:00:00 \n"
                + "#SBATCH -p msismall\n")
        script.write("#SBATCH --nodes=1 \n"
                + "#SBATCH --ntasks=1 \n"
                + "#SBATCH --cpus-per-task=1 \n"
                + "#SBATCH --mem=60gb \n"
                + "#SBATCH --mail-type=ALL \n"
                + "#SBATCH --mail-user=marlo072@umn.edu\n"
                + "#SBATCH --job-name GeneratingSeed" + seed + "\n"
                + "#SBATCH -o /scratch.global/marlo072MCCU/SimulationLOGS/%j.Seed"+seed + ".out\n"
                + "#SBATCH -e /scratch.global/marlo072MCCU/SimulationLOGS/%j.Seed" +seed + ".err\n"
                + "\n"
                + "chmod ug+wrx /scratch.global/marlo072MCCU/SimulationLOGS/* \n\n"
                + "conda activate Python2\n"
                + "newgrp durwa004\n\n"
                + "cd /panfs/jay/groups/6/durwa004/marlo072/neat-genreads \n\n"
                + "python genReads.py \\\n"
                + "\t-r /panfs/jay/groups/6/durwa004/shared/FILESFORSIMULATIONS/goldenPath.Ec_build-3.0_wMSY.fa \\\n"
                + "\t-R " +frag +" \\\n"
                + "\t-o /scratch.global/marlo072/Simulating90/Genomes/Seed" +seed+ " \\\n"
                + "\t-c "+coverage+ " \\\n"
                + "\t-e /scratch.global/marlo072/Simulating90/Models/ErrorModels/" +horse+"ErrorModel \\\n"
                + "\t-M 0 \\\n"
                + "\t-v /scratch.global/marlo072/Simulating90/SeedVCFs/Seed" +seed+"Golden.vcf \\\n"
                + "\t--pe-model /scratch.global/marlo072/Simulating90/Models/FragLengthModels/"+horse+"/fraglen.p \\\n"
                + "\t--gc-model /scratch.global/marlo072/Simulating90/Models/GCModels/"+horse+"GCModel \\\n"
                + "\t--vcf \\\n"
                +"\t--bam \\\n"
                + "\t --gz"
                + "\n\n"
                + "chmod ug+wrx /scratch.global/marlo072/Simulating90/Genomes/*")

##Top bash part (try to put error and output files in my scratch or in my own space so I can watch it)
##Change group
##Sample specific part
