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

### Preview one of the file field categories

```$bioawk -c gff '{print $seqname}' dmel-all-r6-48.gff.gz | head```

### Obtain total incidences of each feature type ordered from most to least common

### First: How many types of features?

```$bioawk -c gff '{types[$feature]++} END {print length(types)}' dmel-allr6.48.gff.gz```

Answer: 58 distinct feature types

### Then proceed to obtain total incidence of each feature type

```$bioawk -c gff '{count[$feature]++} END{ for (f in count) print count[f], f}' dmel-all-r6.48.gff.gz | sort -rn ``` 
Answer:
18133981        match_part

8324113 match

1801236

581034  orthologous_to
254584  oligo

240700  TF_binding_site

208999  polypeptide_region

149517  paralogous_to

133643  RNAi_reagent

114977  CDS

102142  polyA_site

85590   exon

72062   intron

72009   transposable_element_insertion_site

71483   exon_junction

50389   BAC_cloned_genomic_insert

49022   region

43461   TSS

38366   five_prime_UTR

30799   polypeptide

30799   mRNA

28587   three_prime_UTR

25030   regulatory_region

18667   orthologous_region

17896   gene

14095   PCR_product

11413   sgRNA

9791    chromosome_breakpoint

9042    repeat_region

8062    origin_of_replication

7683    insulator

6629    point_mutation

5898    transposable_element

5728    chromosome_band

3053    ncRNA

2870    deletion

2004    modified_RNA_base_feature

1870    golden_path

1831    insertion_site

1401    protein_binding_site

978     syntenic_region

946     delins

855     rescue_region

485     miRNA

365     pseudogene

312     tRNA

300     snoRNA

262     pre_miRNA

246     tandem_repeat

159     insertion

158     MNV

115     rRNA

70      complex_substitution

68      sequence_variant

54      sequence_alteration

32      snRNA

7       mature_protein_region

1       DNA_motif


### Obtain total number of genes per each of the following chromosomes

 X ,Y, 2L 2R, 3L, 3R and 4

```$bioawk -c gff '{print $seqname}' dmel-all-r6.48.gff.gz | sort | uniq -c |sort -rn | head```

Answer:

X: 4725809
Y: 52116
2L: 5254139
2R: 5389966
3L: 5901714
3R:7123617
4: 438309



