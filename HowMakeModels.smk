#This script should be able to both make all the appropriate models and also generate variant VCFs.
import pandas as pd
Horses = open('SimulatingHorsesModels.txt','r')
HorseIDs = Horses.read()
HorseIDs = HorseIDs.split('\n')
Empty = HorseIDs.pop()
Seeds = open('Seedlist.txt','r')
SeedList = Seeds.read()
SeedList = SeedList.split('\n')
Empty = SeedList.pop()
rule all:
	input:
		expand(
			"/scratch.global/marlo072/Simulating90/Models/GCModels/{horse}GCModel",
			horse = HorseIDs
		),
		expand(
			"/scratch.global/marlo072/Simulating90/Models/FragLengthModels/{horse}/{horse}.done",
			horse = HorseIDs
		),
		expand(
			"/scratch.global/marlo072/Simulating90/Models/ErrorModels/{horse}ErrorModel",
			horse = HorseIDs
		),
		expand(
			"/scratch.global/marlo072/Simulating90/BCFStats/Seed{Seed}InsertionStats.txt",
			Seed=SeedList
		)	



include: "variantvcfs.smk"
include: "makingmodels.smk"
