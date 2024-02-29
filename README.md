Dependencies needed to run code:
from Biopython use Entrez and SeqIO to search NCBI database and retrieve sequences
os to pass commands from python to the terminal
numpy to analyze data
Kallisto for sequence quanitification
R studio for running Sleuth program
Sleuth to compareexperimental conditions
Glob to idenitfy files and retrieve their locatrion
Pandas for data analysis
from pathlib use Path to identify files and retrieve their location

Step 1: Choose whether to analyze sample data or all data. If choosing sample data, download the sample data files from the previous page.
If choosing all data, follow the steps below.
Download the SRA files using the wget command followed by the links below
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
For example, downloading the first SRA file: wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

After all the files have been downloaded, use "fastq-dump -I --split-files" to split the files into paired end fastq files
For example, converting the first SRA file: fastq-dump -I --split-files SRR5660030

Following track 1, we will perform differential expression.
Before running the code, be sure the "Sleuth.R" file is in your current working directory.
