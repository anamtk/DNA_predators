# Load biokit
module load biokit
# Load biopython
module load biopython-env

# Concatenate all fasta files into one
cat Greenland_* > Spider2_data.fas

# Fix fasta headers: remove white space, add barcodelabel tag
sed 's/ /_/g' Spider2_data.fas > fix1_Spider2_data.fas
sed 's/mytag=/;barcodelabel=/g' fix1_Spider2_data.fas > fix2_Spider2_data.fas
sed 's/|/;/g' fix2_Spider2_data.fas > fix3_Spider2_data.fasta

# Trim for forward primer until trimmed 0 times
# Cutadapt can be executed as one-liner using option: -n X (X is number of times repeated)
cutadapt -g LCO=CAACAAATCATAAAGATATTGG --minimum-length 100 -e 0.2 fix3_Spider2_data.fasta > fix3_Spider2_data_trim1.fasta
cutadapt -g LCO=CAACAAATCATAAAGATATTGG --minimum-length 100 -e 0.2 fix3_Spider2_data_trim1.fasta > fix3_Spider2_data_trim2.fasta
cutadapt -g LCO=CAACAAATCATAAAGATATTGG --minimum-length 100 -e 0.2 fix3_Spider2_data_trim2.fasta > fix3_Spider2_data_trim3.fasta

# Trim for reverse primer until trimmed 0 times
cutadapt -b MLepR1_rev=GAAAATGGAGCTGGAAC --minimum-length 100 -e 0.2 fix3_Spider2_data_trim3.fasta > fix3_Spider2_data_trim4.fasta
cutadapt -b MLepR1_rev=GAAAATGGAGCTGGAAC --minimum-length 100 -e 0.2 fix3_Spider2_data_trim4.fasta > fix3_Spider2_data_trim5.fasta
cutadapt -b MLepR1_rev=GAAAATGGAGCTGGAAC --minimum-length 100 -e 0.2 fix3_Spider2_data_trim5.fasta > fix3_Spider2_data_trim6.fasta
cutadapt -b MLepR1_rev=GAAAATGGAGCTGGAAC --minimum-length 100 -e 0.2 fix3_Spider2_data_trim6.fasta > fix3_Spider2_data_trim7.fasta
cutadapt -b MLepR1_rev=GAAAATGGAGCTGGAAC --minimum-length 100 -e 0.2 fix3_Spider2_data_trim7.fasta > fix3_Spider2_data_trim8.fasta

# Count sequence lengths
fastalength fix3_Spider2_data_trim8.fasta | cut -d " " -f 1 | sort | uniq -c

# Discard shorter than 300 and truncate at 300
fasta_select_len 300 fix3_Spider2_data_trim8.fasta > min300_reads.fasta
fastx_trimmer -l 300 -i min300_reads.fasta -o global_300bp.fasta

# Dereplication at full length similarity
usearch -derep_fulllength global_300bp.fasta -output derep_global_300bp.fasta -sizeout -uc global_300bp_14Oct2014.uc

# Sorting by haplotype size (and optionally discarding uniques using option: -minsize 2) 
usearch -sortbysize derep_global_300bp.fasta -output sortabundance_derep_global_300bp.fa

# TESTING DATA
usearch -sortbysize derep_global_300bp.fasta -output test1.fa
grep "^>" test1.fa | wc -l
usearch -sortbysize derep_global_300bp.fasta -output test2.fa -minsize 2
grep "^>" test2.fa | wc -l
usearch -cluster_otus test2.fa -otus otus_test2.fa -fulldp
python /wrk/ekrates/drive5_py/fasta_number.py otus_test2.fa OTU_ > label_otus_test2.fa
# TEST ENDS

# Clustering to OTUs at 97% threshold (= default) using full dynamic programming to find alignments with with the maximum possible score
usearch -cluster_otus sortabundance_derep_global_300bp.fa -otus otus.fa -fulldp

# Renaming OTUs
python /wrk/ekrates/drive5_py/fasta_number.py otus.fa OTU_ > label_otus.fa

# Mapping OTUs using global search
usearch -usearch_global sortabundance_derep_global_300bp.fa -db label_otus.fa -strand plus -id 0.97 -uc readmap.uc

# Building a table from OTU mappings
python /wrk/ekrates/drive5_py/uc2otutab.py readmap.uc > readmap.tab

# Blast OTUs to local reference database
# Note: pb blast creates database from a fasta file, that is, no need to build one
pb blastn -dbnuc fixed2_Greenland_Ref_Database_25Aug2014.fas -query label_otus.fa -out label_otus_outfmt6.txt -outfmt 6

# Retrieving species identifications from BOLD using BOLD retriever
nohup python /homeappl/home/ekrates/appl_taito/bold_retriever/bold_retriever/bold_retriever.py -f label_otus.fa -db COX1