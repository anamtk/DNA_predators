#October 21, 2019
# K. Seltmann
# takes a list of names and returns the number of records for a specific gene region based on the plant list and downloads all of the associated fasta files
# requires Cal-plants.txt
# to change the gene, change gene region mentioned in for loop (line31), downloads fasta files based on a list of plants and gene region

#import Entrez
from Bio import Entrez

#change to be your email. NCBI API needs this to work
Entrez.email = "seltmann@ccber.ucsb.edu"

# countGene function that creates a text file of the number of records on ncbii for the plants we are interested in
def countGene(name, gene):
    handle = Entrez.esearch(db='nucleotide', term = [name + "[Orgn] AND " + gene + "[Gene]"])
    record = Entrez.read(handle)
    print(record)
    return record

#function to return fasta file based on single ID
def singleEntry(singleID):   #the singleID is the accession number
    handle = Entrez.efetch(db='nucleotide',id=singleID, rettype = 'fasta', retmode= 'text')
    f = open('%s.fasta' % singleID, 'w')
    f.write(handle.read())
    handle.close()
    f.close()
	
# get list of plant names and put it in a list
text_file = open("Cal-plants.txt", "r")
lines = text_file.readlines()
print(lines)
text_file.close()

# open a file to write to
f = open("plantCounts-downloadFasta.txt", "w")

#go through list and pass to countGene function
for plantNames in lines:
    name = plantNames
    gene = 'rbcl' #change gene region here
    record = countGene(name,gene)
    genes = record["IdList"]
    for g in genes:
        singleEntry(g)
    print(genes)
    recordRow = (name.strip(),",", gene,",",record["Count"]+'\n')
    #print(record)
    f.writelines(recordRow)
	
#close file
f.close()