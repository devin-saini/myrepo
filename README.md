Download the SRA files using the wget command followed by the links below
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
For example, downloading the first SRA file: wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

After all the files have been downloaded, use "fastq-dump -I --split-files" to split the files into paired end fastq files
For example, converting the first SRA: fastq-dump -I --split-files SRR5660030

Following track 2, use Bowtie 2 to create an index for HCMV. First download the genome using:
datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report
After the file has downloaded, unzip the file using: unzip ncbi_dataset.zip
Build index using: bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV

Write out reads that map for D1 2DPI:
bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV-D1-2DPI-map.sam --al-conc-gz SRR5660030_mapped_%.fq.gz
