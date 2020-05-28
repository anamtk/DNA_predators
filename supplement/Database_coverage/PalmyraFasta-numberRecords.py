#October 21, 2019
# from K. Seltmann "PlantFasta-numberRecords.py
# edited by A. Miller-ter Kuile November 2019
# takes a list of names and returns the number of records for a specific gene region based on the species list
# requires pal_species_nov19.txt
# to change the gene, change gene region mentioned in for loop (line31)

#import Entrez
from Bio import Entrez

#change to be your email. NCBI API needs this to work
Entrez.email = "ana00@ucsb.edu"

# countGene function that creates a text file of the number of records on ncbii for the plants we are interested in
def countGene(name, gene):
    handle = Entrez.esearch(db='nucleotide', term = [name + "[Orgn] AND " + gene + "[All Fields]"])
    #unsure about if should be searching specifically by Gene or by All Fields
    #  handle = Entrez.esearch(db='nucleotide', term = [name + "[Orgn] AND " + gene + "[Gene]"])
    record = Entrez.read(handle)
    print(record)
    return record
	
# get list of plant names and put it in a list
text_file = open("pal_species_nov19.txt", "r") #this is the file with organism names in a list
lines = text_file.readlines()
print(lines)
text_file.close()

# open a file to write to
#change name of file based on each gene region
f = open("counts_pal_[gene]_nov19.txt", "w") #this is the file you are creating that will have your names and genes with numbers

#go through list and pass to countGene function 
for plantNames in lines:
    name = plantNames
    gene = 'coi' #change gene region here, did both co1 and coi with different outputs
    record = countGene(name,gene)
    recordRow = (name.strip(),",", gene,",",record["Count"]+'\n')
    print(recordRow)
    f.writelines(recordRow)
	
#close file
f.close()
