# From https://www.gencodegenes.org/human/; Release 47 (GRCh38.p14); selected based on guidance in star documentation page 6 (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

# Download the human reference transcriptome (FASTA)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# Download the GENCODE annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz

# Unzip the downloaded files
gunzip -c GRCh38.primary_assembly.genome.fa.gz > reference/GRCh38.primary_assembly.genome.fa
gunzip -c gencode.v47.primary_assembly.annotation.gtf.gz > reference/gencode.v47.primary_assembly.annotation.gtf
