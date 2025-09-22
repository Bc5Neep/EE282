# Homework #3

## Initiate run on dedicated cluster

srun -A class_ee282 --pty --x11 bash -i

# Initialize mamba
mamba activate ee282

# Retreive June 2022 release of  D.melanogaster
wget https://ftp.flybase.net/releases/FB2022_05/dmel_all_CDS_r6.48/fasta/

### Unpack file and inspect elements 
faSize  dmel-all-CDS-r6.48.fasta.gz
# View tab delimited file
faSize -tab dmel-all-CDS-r6.48.fasta.gz

# Check file integrity
md5sum dmel-all-CDS-r6.48.fasta.gz

# Calculate summaries of genome

# Summarize Annotation
wget https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/

# Check file Integrity
md5sum dmel-all-r6.48.gff.gz

# Preview file categories

bioawk -c help

# Preview one of the file field categories

bioawk -c gff '{print $seqname}' dmel-all-r6.48.gff.gz | head 


# Obtain total number of features of each type ordered from most to least common
# First: How many features?

bioawk -c gff '{ types[$feature]++ } END {print length (types)}' dmel-all-r6.48.gff.gz

#Then total number of features of each type
bioawk -c gff '{count[$feature]++} END {for (f in count) print count[f],f}' dmel-all-r6.48.gff.gz | sort -rn  

# Obtain total number of genes per each following chromosome 
bioawk -c gff '{print $seqname}' dmel-all-r6.48.gff.gz | sort | uniq -c |sort -rn | head



