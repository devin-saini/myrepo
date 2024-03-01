Dependencies needed to run code:
from Biopython import Entrez and SeqIO to analyze NCBI files, search NCBI database, and retrieve sequences
os to pass commands from python to the terminal
numpy to analyze data
Kallisto for sequence quanitification
R studio for running Sleuth program
Sleuth to compare experimental conditions
Glob to idenitfy files and retrieve their location
Pandas for data analysis
from pathlib import Path to identify files and retrieve their location

Ensure that all steps below are performed in the same directory.

Step 1: Choose whether to analyze sample data or all data. If choosing sample data, download the sample data files from the previous page and begin at Step 2.
If choosing all data, in your current working directory, follow the steps below.
Download the SRA files using the wget command followed by the links below
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
For example, downloading the first SRA file: wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

After all the files have been downloaded, use "fastq-dump -I --split-files" to split the files into paired end fastq files
For example, converting the first SRA file: fastq-dump -I --split-files SRR5660030

Step 2: Following track 1, we will perform differential expression.
Before running the code, be sure the "Sleuth.R" file and ALL FASTQ files are in the same directory as the code.

After ensuring all appropriate files are in the same directory and all fastq files end with ".fastq", run the code by entering "python Python_Pipeline.py" in the command line.
