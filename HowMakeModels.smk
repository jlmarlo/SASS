#This script should be able to both make all the appropriate models and also generate variant VCFs.
with open(config['samples'],'r') as f:
    HorseIDs = [line.strip() for line in f if line.strip()]


rule all:
    input:
        expand(
            "outputs/Models/GCModels/{sample}GCModel",
            sample = HorseIDs),
        #expand(
        #    "outputs/Models/FragLengthModels/{sample}/{sample}.done",
        #    sample = HorseIDs),
        #expand(
        #    "outputs/Models/ErrorModels/{sample}ErrorModel",
        #    sample = HorseIDs),
        #"outputs/BCFStats/all.done"


include:"rules/makingmodels.smk"
include: "rules/variantvcfs.smk"
