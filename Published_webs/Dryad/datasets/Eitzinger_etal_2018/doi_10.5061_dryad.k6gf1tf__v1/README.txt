This file explains all variables in the in the datasets that accompany

Assessing changes in arthropod predator-prey interactions through DNA-based gut content analysis – variable environment, stable diet
submitted to Molecular Ecology

Bernhard Eitzinger1, Nerea Abrego1, Dominique Gravel2, Tea Huotari1, Eero J Vesterinen1,3 and Tomas Roslin1,4

1 Department of Agricultural Sciences, Latokartanonkaari 5, FI-00014 University of Helsinki, Finland 
2 Département de biologie, Université de Sherbrooke, 2500 Boulevard de l’Université, Sherbrooke, Quebec, Canada J1K 2R1
3 Biodiversity Unit, University of Turku, Vesilinnantie 5, FI-20014 Turku, Finland
4 Department of Ecology, Swedish University of Agricultural Sciences, Box 7044, 750 07 Uppsala, Sweden

contact: Bernhard Eitzinger nanuuto@gmail.com


###Arthropod community data matrices

Study design 
	hierarchical: sampling unit (plots) replicated temporally 3 times (t2,t3,t4 field data) and 4 times (t1,t2,t3,t4 spider gut content data) respectively
	plot names: 18 plots A,B,C,D,DX,E,F,G,H,I,IX,J,K,L,LX,P,R,S
	plotID: combination of plot names and sampling date

A_1 matrix (animal trap data)	
	abundance data of prey taxa (family level): Individuals/trapday	
	three sampling dates per plot: t2,t3,t4 corresponding to sampling dates in other matrices	
	48 taxa x 54 sampling units
	
A_2 matrix (animal trap data)	
	note: data in A_2 matrix is just re-formated data from A_1 matrix	
	presence/absence data of prey taxa (family level)	
	three sampling dates per plot: t2,t3,t4 corresponding to sampling dates in other matrices	
	48 taxa x 54 sampling units	

B_1 matrix (spider gut content data and predator traits)	
	spider gut content (individual wolf spiders tested for prey taxa, plot level): 47 taxa x 670 sampling units (individual spiders; pardosaID)	
	four sampling dates per plot: t1, t2,t3,t4 corresponding to sampling dates in other matrices	
	pred.size (predator body size; mm)	
	pred.mass(mg freshweight; calculated from from non-linear length-mass regression model)	
	pred.devstage (predator development stage; Immature/Mature)	
	pred.sex (predatorsex; Female/Male)
	
B_2 matrix (spider gut content data)	
	note: data in B_2 matrix is just re-formated data from B_1 matrix	
	spider gut content (mean number of wolf spiders tested positive for prey taxa, plot level): 47 taxa x 72 sampling units	
	four sampling dates per plot: t1, t2,t3,t4 corresponding to sampling dates in other matrices

C matrix (prey traits)	
	47 taxa x 3 traits
	Feeding mode (herbivore, detritivore,carnivore, omnivore,parasitoid)
	Habitat (aboveground, belowground)
	Body mass (g dryweight)

D matrix (environmental data)	
	72 sampling units x 9 covariates
	GPS coordinates (grid zone, UTM_E, UTM_N)
	Time of sampling (t1,t2,t3,t4)
	Elevation (m a.s.l.)
	TempAir (mean air temperature from time tx-1 to tx; degree Celsius; measured at plot level)
	TempSoil (mean soil temperature from time tx-1 to tx; degree Celsius; measured at plot level)
	HumidRel (mean relative humidity from time tx-1 to tx; in %; measured at plot level)
	Veg.mass (dry mass vegetation+litter in g from a soil core with 4 cm diameter; measured at plot level)
	VegCover.dryas (% Dryas plants on 50 x 50 cm plot; measured at plot level)
	Soil.water (% water in 4 cm diameter soil core; measured at plot level)


###Sequence and Bioinformatic datasets

ARC_ZOTUS_Sequences and Zeale_ZOTUS_Sequences
	Contains separate sequences for each primer set clustered to zero-radius OTUs (ZOTUs) using USEARCH algorithm.

ARC_readmap and Zeale_readmap linking trimmed reads and samples
	These readmaps shows read count for each of the two PCR replicates per predator sample. The readmap is constructed using USEARCH software. Names of the 		samples correspond to data in the matrix B_1

ARC and Zeale final readmap
	These readmaps contains the final read counts for each predator sample.

Command pipeline used for bioinformatics
The file contains all the commands used in this study to generate the final data. The bioinformatics was carried out at servers on CSC - IT Center for Science in Finland <www.csc.fi>.

	



