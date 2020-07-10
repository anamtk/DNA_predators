# Merge paired-end reads, while relabelling each reads with sample name and #,
# then filtering for quality and discarding short reads
usearch -fastq_mergepairs *_R1_001.fastq -fastq_maxdiffs 10 -fastqout merged.fq \
-fastaout_notmerged_fwd merged_fwd.notmerged \
-fastaout_notmerged_rev merged_rev.notmerged -relabel @
# Trim for quality using max_ee-rate 1 for the whole read
usearch -fastq_filter merged.fq -fastq_maxee_rate 1 -fastaout reads.fasta
# Remove PCR primers after merging and add trimmed primer information
# in the end of the read header
# Anchor forward primers to 5' end using ^ and reverse primers to 3' end using $
# NOTE: When removing primers AFTER merging, the reverse primers need to 
# be reverse complements!
# Cut forward primers from 5' end and reverse primers from the 3' at the same time
# Prepare a file named 'primers.fa' containing all primer pairs:
# These are linked primers and 5' primer is anchored by default
>ZBJ-ArtF1c_ZBJ-ArtR2c
AGATATTGGAACWTTATATTTTATTTTTGG...GGAGGATTTGGWAATTGATTAGTW
>ARC-F3_ARC-R6 
GCTTTYCCYCGAATAAATAATWTAAG...AGGWTGAACWGTTTAYCCTCC
# Trim primer pairs using the file 'primers.txt' as primer input source
# Too short reads with '--too-short-output' option into a file
# Length limit for Zeale=100 and for ARC=50
cutadapt -a file:primers.fa --minimum-length=50 --error-rate=0.2 \
--untrimmed-output=untrimmed.fa --too-short-output tooShort.fa \
-o labeled_{name}.fasta reads.fasta 
# Dereplicate and relabel 
usearch -fastx_uniques labeled_ARC-F3_ARC-R6.fasta -sizeout -relabel Uniq \
-fastaout ARC_uniq.fa
usearch -fastx_uniques labeled_ZBJ-ArtF1c_ZBJ-ArtR2c.fasta -sizeout -relabel Uniq \
-fastaout Zeale_uniq.fa
# Chimera filtering using UNOISE either do zotus or otus
usearch -unoise3 ARC_uniq.fa -zotus ARC_zotus.fa
usearch -unoise3 Zeale_uniq.fa -zotus Zeale_zotus.fa
# Rename ZOTUs as OTUs for compatibility
sed -i 's/Zotu/Otu/g' ARC_zotus.fa
sed -i 's/Zotu/Otu/g' Zeale_zotus.fa
# Do 97% OTUs and remove chimerics
usearch -cluster_otus ARC_uniq.fa -otus ARC_otus.fa -relabel Otu
usearch -cluster_otus Zeale_uniq.fa -otus Zeale_otus.fa -relabel Otu
# Build (z)otutabs
usearch -otutab labeled_ARC-F3_ARC-R6.fasta -sample_delim . -otus ARC_zotus.fa \
-otutabout ARC_zotus.tab
usearch -otutab labeled_ZBJ-ArtF1c_ZBJ-ArtR2c.fasta -sample_delim . -otus Zeale_zotus.fa \
-otutabout Zeale_zotus.tab
usearch -otutab labeled_ARC-F3_ARC-R6.fasta -sample_delim . -otus ARC_otus.fa \
-otutabout ARC_otus.tab
usearch -otutab labeled_ZBJ-ArtF1c_ZBJ-ArtR2c.fasta -sample_delim . -otus Zeale_otus.fa \
-otutabout Zeale_otus.tab
# Use bold_retriever to get bold identifications for otus
python bold_retriever.py -f ARC_zotus.fa -db COX1
python bold_retriever.py -f Zeale_zotus.fa -db COX1

# BLAST against own database with output format 6
blastn -db Greenland_BLASTdb -query ARC_zotus.fa -out ARC_Greenland_BLAST.txt \
-outfmt "6 qseqid stitle sblastnames sseqid pident length mismatch gapopen qstart \
qend sstart send evalue bitscore staxids sscinames sskingdoms"
# BLAST to GenBank nucleotide database with output format 6
blastn -db nt -query ARC_zotus.fa -out ARC_Greenland_pbBLAST_custom.txt \
-outfmt "6 qseqid stitle sblastnames sseqid pident length mismatch gapopen qstart \
qend sstart send evalue bitscore staxids sscinames sskingdoms"
# Fields (comma-separated):
# Query seq-id, Subject title, Subject blast Name(s), Subject seq-id,
# Percentage of identical matches, Alignment length, Number of mismatches,
# Number of gap openings, Start of alignment in query, End of alignment in query,
# Start of alignment in subject, End of alignment in subject, Expect value, Bit score,
# Subject tax ids, Subject sci names, subject super kingdoms