from Bio import Entrez
from Bio import SeqIO
import os
import numpy as np
import glob
import pandas as pd
from pathlib import Path

#identify all fastq files and sort them
fastq = glob.glob('*.fastq')
fastq = sorted(fastq)

#make file paths for all fastq files
fastq_filepath = []
for f in fastq:
  filepath = Path(f).absolute()
  fastq_filepath.append(filepath)

os.mkdir("PipelineProject_Devinjeet_Saini") #make directory to store all files made

#move sleuth file to pipeline  directory
sleuth_r_path = Path('Sleuth.R').absolute() #path to sleuth code
new_spot = Path("PipelineProject_Devinjeet_Saini").absolute() #path to pipeline directory
new_spot = f'{new_spot}/Sleuth.R' #path where to store sleuth code
os.rename(sleuth_r_path,new_spot) #move sleuth code

os.chdir("PipelineProject_Devinjeet_Saini") #move into directory

#get genbank file from of HCMV
accession = 'NC_006273.2'
Entrez.email = 'dsaini@luc.edu'
handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb',retmode = 'text')
record = SeqIO.read(handle,'genbank')

#count coding sequences in genbank file
cds_count = 0
for feature in record.features:
    if feature.type == "CDS":
       cds_count += 1

#make new file and add number of coding sequences in HCMV
with open ('PipelineProject.log','w') as f:
  f.write(f'The HCMV genome (NC_006273.2) has {cds_count} CDS.\n')

#make fasta file of HCMV cds
with open ('HCMV_cds.fasta','w') as outfile:
   for feature in record.features:
      if feature.type == "CDS":
         protein_id = feature.qualifiers.get('protein_id')
         start = feature.location.start #identify start position of CDS
         end = feature.location.end #identify end position of CDS
         sequence = record.seq[start:end] #identify CDS sequence
         outfile.write(f">{protein_id[0]}\n{sequence}\n")

hcmv = "HCMV.idx" #where to store index
hcmv_cds = Path('HCMV_cds.fasta').absolute() #path to reference trancriptome
os.system(f'kallisto index -i {hcmv} {hcmv_cds}') #make kallisto of HCMV CDS

index = Path('HCMV.idx').absolute() #path to reference trancriptome

tpm_output = ['D1-2DPI','D1-6DPI','D3-2DPI','D3-6DPI'] #Make list of all conditions

current_dir = os.getcwd() #identify current working directory

#Find tpm of values of samples
output = tpm_output[0]
output_path = os.path.join(current_dir,output) #file path of output
read1 = fastq_filepath[0] #file path of paired end read 1
read2 = fastq_filepath[1] #file path of paired end read 2
os.system(f'time kallisto quant -i {index} -o {output} -b 30 -t 4 {read1} {read2}')

output = tpm_output[1]
output_path = os.path.join(current_dir,output) #file path of output
read1 = fastq_filepath[2] #file path of paired end read 1
read2 = fastq_filepath[3] #file path of paired end read 2
os.system(f'time kallisto quant -i {index} -o {output} -b 30 -t 4 {read1} {read2}')

output = tpm_output[2]
output_path = os.path.join(current_dir,output) #file path of output
read1 = fastq_filepath[4] #file path of paired end read 1
read2 = fastq_filepath[5] #file path of paired end read 2
os.system(f'time kallisto quant -i {index} -o {output} -b 30 -t 4 {read1} {read2}')

output = tpm_output[3]
output_path = os.path.join(current_dir,output) #file path of output
read1 = fastq_filepath[6] #file path of paired end read 1
read2 = fastq_filepath[7] #file path of paired end read 2
os.system(f'time kallisto quant -i {index} -o {output} -b 30 -t 4 {read1} {read2}')

#add headers for TPM
with open ('PipelineProject.log','a') as f:
  f.write('\n')
  f.write('sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n')

#Find all kallisto files
globby = glob.glob('D*-*DPI')
globby = sorted (globby)

#make list of all samples and conditions
sample = ['SRR566030','SRR566033','SRR566044','SRR566045']
condition = ['2dpi','6dpi','2dpi','6dpi']

#for each kallisto file, identify desired tpm values and add those values to log file
for d in range (len(globby)):
  os.chdir(globby[d]) #move into kallisto file
  data = np.loadtxt('abundance.tsv', delimiter="\t", skiprows=1, usecols=(-1,)) #identify tpm column
  min_tpm = np.min(data).round(2) #calculate minimum
  median_tpm = np.median(data).round(2) #calculate median
  mean_tpm = np.mean(data).round(2) #calculate mean
  max_tpm = np.max(data).round(2) #calculate maximum
  os.chdir('..') #move back 1 directory
  with open ('PipelineProject.log','a') as f: #open log file and write to it
    f.write(f"{sample[d]}\t{condition[d]}\t{min_tpm}\t{median_tpm}\t{mean_tpm}\t{max_tpm}\n")

#open log file and add headers row for differentiation values
with open ('PipelineProject.log','a') as f: 
  f.write('\n')
  f.write('target_id\ttest_stat\tpval\tqval\n')

#make table for sleuth code
with open('sleuth_table.txt','w') as f:
  f.write ('sample condition path\n')
  f.write('SRR5660030 2dpi D1-2DPI\n')
  f.write('SRR5660044 2dpi D3-2DPI\n')
  f.write('SRR5660033 6dpi D1-6DPI\n')
  f.write('SRR5660045 6dpi D3-6DPI')

#make emtpy lists to hold values
target_id = []
test_stat = []
pval = []
qval = []

os.system('Rscript Sleuth.R') #Run sleuth on kallisto files

#open sleuth output and identify each value for each transcript
with open ('sleuth_results.txt','r') as f:
  next(f)
  for line in f:
    values = line.strip().split()
    target_id.append(values[0])
    test_stat.append(values[3])
    pval.append(values[1])
    qval.append(values[2])

#open log file and write details for each row
with open ('PipelineProject.log','a') as f:
  for i in range (len(target_id)):
    f.write(f"{target_id[i]}\t{test_stat[i]}\t{pval[i]}\t{qval[i]}\n")

cds_id = target_id[0] #identify most differentially expressed protein id

#get fasta file of most differentially expressed protein
Entrez.email = 'dsaini@luc.edu'
handle = Entrez.efetch(db='protein', id=cds_id, rettype='fasta')
protein_fasta = SeqIO.read(handle,'fasta')

#download fasta file of most differentiall expressed protein
with open('protein.fasta','w') as f:
  SeqIO.write(protein_fasta,f,'fasta')

#download all nucleotides for Betaherpesvirinae family
term = 'Betaherpesvirinae'
os.system(f'datasets download virus genome taxon {term} --include genome')
os.system('unzip ncbi_dataset.zip')

#make local nucleotide database for Betaherpesvirinae
nuc_fasta = Path('ncbi_dataset/data/genomic.fna').absolute()
os.system(f'makeblastdb -in {nuc_fasta} -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl')

#blast protein against Betaherpesvirinae nucleotides
blast_input = 'protein.fasta'
output_file = 'myresults.csv'
os.system(f'tblastn -query {blast_input} -db Betaherpesvirinae -out {output_file} -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"')

#create header values for blast output
with open ('PipelineProject.log','a') as f:
  f.write('\n')
  f.write('sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')

#make dataframe of blast
df = pd.read_csv('myresults.csv', nrows=10, usecols=range(10),header=None)

#add dataframe to log file
df.to_csv('PipelineProject.log', mode='a', sep='\t', index=False, header=False)
