# Homework #3

## Initiate run on dedicated cluster

```$srun -A class_ee282 --pty --x11 bash -i```

## Initialize mamba

```mamba activate ee282```

### Retreive June 2022 release of [*Drosophila Melanogaster*](https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-CDS-r6.48.fasta.gz) fasta file from Flybase

```wget https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/```

### Unpack file and inspect elements

```faSize https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/```

### View tab delimited file

```faSize -tab https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/```

### Check file integrity

```md5sum --ignore missing -c md5sum.txt```

### Calculate summaries of genome

1. Total number of nucleotides Answer: 61223799
2. Total number of Ns Answer: 0
3. Total number of sequences Answer: 30799

### Summarize Annotation

Obtain [.gtf file of *Drosophila Melanogaster*](https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/)
from Flybase

Check file Integrity

```md5sum dmel-all-r6.48.gff.gz```

### Preview file categories

$bioawk -c help

Preview one of the file field categories

```$bioawk -c gff '{print $seqname}' dmel-all-r6-48.gff.gz | sort -ru | head```

### Obtain top 10 list of features that appear most frequently in each category

``` $bioawk -c gff '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' dmel-all-r6.48.gff.gz | sort | uniq -c | sort -rn | head ```

### Obtain total number of genes per chromosome of the following

 X ,Y, 2L 2R, 3L, 3R and 4

```$bioawk -c gff '{print $seqname}' dmel-all-r6.48.gff.gz | sort | uniq -c |sort -rn | head```
