# convert bam to fastq
./bamtofastq-1.2.0 possorted_genome_bam.bam ./exp1/
# start processing fastq
conda create --name kb_env python=3.7
pip install kb-python
conda install -c bioconda bustools

# download ref
kb ref -d linnarsson -i index.idx -g t2g.txt -c1 spliced_t2c.txt -c2 unspliced_t2c.txt

# kb count --h5ad -i index.idx -g t2g.txt -x 10xv2 -o MCF10Atimecourse \
# -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 8 \
# MCF10Atimecourse_S1_L001_R1_001.fastq.gz.1 \
# MCF10Atimecourse_S1_L001_R2_001.fastq.gz.1

# direct download google file
# wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=11wWUBqjEmy-6byjv6dbrLDtmVL7tcPkD' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=11wWUBqjEmy-6byjv6dbrLDtmVL7tcPkD" -O exp1.bam && rm -rf /tmp/cookies.txt
