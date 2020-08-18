# convert bam to fastq
./bamtofastq_linux possorted_genome_bam.bam ./exp1/
# start processing fastq
conda create --name kb_env python=3.7
pip install kb-python
conda install -c bioconda bustools

# download ref
# kb ref -d linnarsson -i index.idx -g t2g.txt -c1 spliced_t2c.txt -c2 unspliced_t2c.txt
kb ref -d human -i index.idx -g t2g.txt -c1 spliced_t2c.txt -c2 unspliced_t2c.txt
# kb count --h5ad -i index.idx -g t2g.txt -x 10xv2 -o MCF10Atimecourse \
# -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 8 \
# MCF10Atimecourse_S1_L001_R1_001.fastq.gz.1 \
# MCF10Atimecourse_S1_L001_R2_001.fastq.gz.1
kb count --verbose --h5ad -m 16G -i index.idx -g t2g.txt -x 10xv2 -o kazu_mcf10a_10xv2_linRef_lamanno -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 16\
 ./exp1/RXhi10056_MissingLibrary_1_H3T27BCX2/*fastq.gz\
 ./exp1/RXhi10056_MissingLibrary_1_HCKK5BCX2/*fastq.gz

kb count --verbose --h5ad -m 16G -i index.idx -g t2g.txt -x 10xv2 -o kazu_mcf10a_10xv2_samtools_linRef -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 16\
kazu_exp1.fastq.gz
# direct download google file
# wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=11wWUBqjEmy-6byjv6dbrLDtmVL7tcPkD' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=11wWUBqjEmy-6byjv6dbrLDtmVL7tcPkD" -O exp1.bam && rm -rf /tmp/cookies.txt

kb count --verbose --h5ad -m 8G -w kazu_whitelist_barcodes.txt -i index.idx -g t2g.txt -x 10xv3 -o kazu_mcf10a_10xv3_kazuWhitelist -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 8\
 ./exp1/RXhi10056_MissingLibrary_1_H3T27BCX2/*fastq.gz\
 ./exp1/RXhi10056_MissingLibrary_1_HCKK5BCX2/*fastq.gz

samtools fasta ../possorted_genome_bam.bam > kazu_exp1.fasta
gzip kazu_exp1.fastq

kb ref -i kazu_index.idx -g kazu_t2g -f1 kazu_cdna_fasta -f2 kazu_intron_fasta -c1 kazu_c1 -c2 kazu_c2 ../possorted_genome_bam.bam


bustools correct -w 10xv2_whitelist.txt output.bus -o corrected_output.bus
bustools text corrected_output.bus -o corrected_output.txt




# velocyto pipeline
# installation
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install velocyto

