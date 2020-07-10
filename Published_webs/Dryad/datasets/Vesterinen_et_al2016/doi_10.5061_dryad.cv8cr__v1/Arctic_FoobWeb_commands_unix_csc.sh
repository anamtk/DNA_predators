# Load biokit
module load biokit
# Load biopython
module load biopython-env
# Load fastqc
module load fastqc

# Strip linker-tags plus barcodes and rename all barcoded reads to include individual barcode and also rename reads
# IMPORTANT: Change the barcodes_birds.txt file headers according to the samples!
nohup python /wrk/ekrates/drive5_py/fastq_strip_barcode_relabel2.py Greenland_9Sep2014.fastq GATACGACGTTGTAAAA barcodes_birds.txt birdprey > labeled_reads_birds.fastq
nohup python /wrk/ekrates/drive5_py/fastq_strip_barcode_relabel2.py  Spider_gut_9Sep2014.fastq GATACGACGTTGTAAAA barcodes_spiders.txt spiderRead > labeled_reads_spiders.fastq

# Discard shorter than 100 using FASTX Toolkit since there may be some blank sequences causing problems later
nohup fastx_clipper -v -l 100 -i labeled_reads_birds.fastq -o shortOut_birds.fastq -Q 33
nohup fastx_clipper -v -l 100 -i labeled_reads_spiders.fastq -o shortOut_spiders.fastq -Q 33

# Splitting reads into forward and reverse reads
nohup perl /wrk/ekrates/split.pl shortOut_birds
# _frw.fastq = forward reads
# _rvr.fastq = reverse reads
# _vsh.fastq = very short reads
# _oth.fastq = other reads

nohup sed 's/barcodelabel=/barcodelabel=frw_/g' shortOut_birds_frw.fastq > frw_reads_birds.fastq
nohup sed 's/barcodelabel=/barcodelabel=rvr_/g' shortOut_birds_rvr.fastq > rvr_reads_birds.fastq

# FORWARD READS birds
# Cutadapt can be executed as one-liner using option: -n X (X is number of times repeated)
#
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 frw_reads_birds.fastq > frw_birds_trim1.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 frw_birds_trim1.fastq > frw_birds_trim2.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 frw_birds_trim2.fastq > frw_birds_trim3.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 frw_birds_trim3.fastq > frw_birds_trim4.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 frw_birds_trim4.fastq > frw_birds_trim5.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 frw_birds_trim5.fastq > frw_birds_trim6.fastq
#############################################################################################################################
#
# REVERSE READS birds
#
nohup cutadapt -g R-primer=NACTAATCAATTNCCAAATCCTCC --minimum-length 100 -e 0.2 rvr_reads_birds.fastq > rvr_birds_trim1.fastq
nohup cutadapt -g R-primer=NACTAATCAATTNCCAAATCCTCC --minimum-length 100 -e 0.2 rvr_birds_trim1.fastq > rvr_birds_trim2.fastq
nohup cutadapt -g R-primer=NACTAATCAATTNCCAAATCCTCC --minimum-length 100 -e 0.2 rvr_birds_trim2.fastq > rvr_birds_trim3.fastq
nohup cutadapt -g R-primer=NACTAATCAATTNCCAAATCCTCC --minimum-length 100 -e 0.2 rvr_birds_trim3.fastq > rvr_birds_trim4.fastq
nohup cutadapt -g R-primer=NACTAATCAATTNCCAAATCCTCC --minimum-length 100 -e 0.2 rvr_birds_trim4.fastq > rvr_birds_trim5.fastq
nohup cutadapt -g R-primer=NACTAATCAATTNCCAAATCCTCC --minimum-length 100 -e 0.2 rvr_birds_trim5.fastq > rvr_birds_trim6.fastq
#############################################################################################################################
#
# FORWARD READS spiders
#
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 shortOut_spiders.fastq > spiders_trim1.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 spiders_trim1.fastq > spiders_trim2.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 spiders_trim2.fastq > spiders_trim3.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 spiders_trim3.fastq > spiders_trim4.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 spiders_trim4.fastq > spiders_trim5.fastq
nohup cutadapt -g F-primer=AGATATTGGAACNTTATATTTTATTTTTGG --minimum-length 100 -e 0.2 spiders_trim5.fastq > spiders_trim6.fastq
#############################################################################################################################

# Chart qualities before starting filtering FastQC
nohup fastqc frw_birds_trim6.fastq
nohup fastqc rvr_birds_trim6.fastq
nohup fastqc spiders_trim6.fastq

# Chart qualities before starting filtering USEARCH
nohup usearch -fastq_chars frw_birds_trim6.fastq -log frw_birds_trim6_chars.log
nohup usearch -fastq_chars rvr_birds_trim6.fastq -log rvr_birds_trim6_chars.log
nohup usearch -fastq_chars spiders_trim6.fastq -log spiders_trim6_chars.log

# Filtering for quality separately for forward and reverse reads
nohup usearch -fastq_filter frw_birds_trim6.fastq -eeout -fastq_maxee 0.6 -fastq_trunclen 160 -fastq_ascii 33 -fastq_qmin 1 -fastq_qmax 45 -fastaout filtered_frw_birds_trim6.fa -fastqout filtered_frw_birds_trim6.fastq
nohup usearch -fastq_filter rvr_birds_trim6.fastq -eeout -fastq_maxee 0.6 -fastq_trunclen 160 -fastq_ascii 33 -fastq_qmin 1 -fastq_qmax 45 -fastaout filtered_rvr_birds_trim6.fa -fastqout filtered_rvr_birds_trim6.fastq
nohup usearch -fastq_filter spiders_trim6.fastq -eeout -fastq_maxee 0.6 -fastq_trunclen 160 -fastq_ascii 33 -fastq_qmin 1 -fastq_qmax 45 -fastaout filtered_spiders_trim6.fa -fastqout filtered_spiders_trim6.fastq

# Chart qualities AFTER quality filtering FastQC
nohup fastqc filtered_frw_birds_trim6.fastq
nohup fastqc filtered_rvr_birds_trim6.fastq
nohup fastqc filtered_spiders_trim6.fastq

# Reverse complement bird reverse reads before concatenating
nohup fastx_reverse_complement -i filtered_rvr_birds_trim6.fa -o revcomp_filtered_rvr_birds_trim6.fa

# Combine all reads into one file since
# This is biologically meaningful if all the reads originate from same system
nohup cat filtered_frw_birds_trim6.fa revcomp_filtered_rvr_birds_trim6.fa filtered_spiders_trim6.fa > all_reads_11Oct2014.fa

# Dereplication at full length similarity
nohup usearch -derep_fulllength all_reads_11Oct2014.fa -output derep_all_reads_11Oct2014.fa -sizeout -uc derep_all_reads_11Oct2014.uc -minuniquesize 2

# Sorting by haplotype size (and optionally discarding uniques using option: -minsize 2) 
nohup usearch -sortbysize derep_all_reads_11Oct2014.fa -minsize 2 -output sortabundance_derep_all_reads_11Oct2014.fa

# Clustering to OTUs at 97% threshold (= default) using full dynamic programming to find alignments with with the maximum possible score
nohup usearch -cluster_otus sortabundance_derep_all_reads_11Oct2014.fa -otus otus.fa -fulldp

# Renaming OTUs
nohup python /wrk/ekrates/drive5_py/fasta_number.py otus.fa OTU_ > label_otus.fa

# MApping OTUs using global search
nohup usearch -usearch_global sortabundance_derep_all_reads_11Oct2014.fa -db label_otus.fa -strand plus -id 0.97 -uc readmap.uc

# Building a table from OTU mappings
nohup python /wrk/ekrates/drive5_py/uc2otutab.py readmap.uc > readmap.tab

#############################################################################################################################
# Build a local reference database using a fasta file: Greenland_Ref_Database_25Aug2014.fas
# module load emboss
# Remove bad characters from fasta and output to fixed_*
# perl /wrk/ekrates/fastafix_removeIUPAC.pl -i Greenland_Ref_Database_25Aug2014.fas
# sed 's/[().]//g' fixed_Greenland_Ref_Database_25Aug2014.fas > fixed2_Greenland_Ref_Database_25Aug2014.fas
#############################################################################################################################
# seqret fixed2_Greenland_Ref_Database_25Aug2014.fas fixed2_Greenland_Ref_Database_25Aug2014_ncbi.fas -osf ncbi
# makeblastdb -in fixed2_Greenland_Ref_Database_25Aug2014_ncbi.fas -out fixed2_Greenland_Ref_Database_25Aug2014_ncbi -parse_seqids -dbtype nucl
# blastdbcmd -db fixed_Greenland_Ref_Database_25Aug2014_ncbi -info
#############################################################################################################################
# Blast OTUs to (newly created) local reference database
# Note: pb blast creates database from a fasta file, that is, no need to build one
# nohup pb blastn -dbnuc fixed2_Greenland_Ref_Database_25Aug2014.fas -query label_otus.fa -out label_otus_outfmt6.txt -outfmt 6
#############################################################################################################################

# Retrieving species identifications from BOLD using BOLD retriever
nohup python /homeappl/home/ekrates/appl_taito/bold_retriever/bold_retriever/bold_retriever.py -f label_otus.fa -db COX1