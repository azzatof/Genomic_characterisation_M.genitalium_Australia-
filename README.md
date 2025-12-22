# Genomic characterisation of *Mycoplasma genitalium* in Victoria, Australia, reveals lineage diversification and drivers of antimicrobial resistance
Francesca Azzato, George Taiaroa , Janath Fernando, Mona L Taouk, Vesna De Petra, Lenka A Vodstrcil, Erica L Plummer, Kerry Raios, Danielle J Ingle, Niamh Meagher, Jacqueline Prestedge, Eloise Williams, Leon Caly, Benjamin P Howden, Shivani Pasricha, Catriona S Bradshaw, Deborah A Williamson

This GitHub repository contains all code used in this study. Sequence reads are available from the NCBI database under BioProjects PRJNA1367946 and PRJEB5172, with accession numbers provided in the Supplementary Dataset accompanying this publication. All analyses performed in this study can be replicated using the code provided, with input names, directory paths, and file names modified as required.


# Quality Control of *Mycoplasma genitalium* sequence reads

### 1. Pre-filtering of Sequences

Local and global sequence have been pre-filtered to select for reads that only mapped to the G37 *Mycoplasma genitalium* G37 reference genome (NC_000908.2) using minimap2 (v2.24) and samtools (v1.10). 

Additionally all sequence reads underwent adapter trimming using Trimmomatic (v 0.39) prior to QC analysis to remove residual sequence adapters.

```
trimmomatic PE -threads 8 -phred33 <sampleID>_R1_fastq.gz <sampleID>_R2_fastq.gz <sampleID>_Trim_R1.fastq.gz tmp1.fastq.gz <sampleID>_Trim_R2.fastq.gz tmp2.fastq.gz ILLUMINACLIP:Adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
### 2.Sequencing depth
Trimmed paried end reads for each genome were aligned to the *M.genitlium* G37 reference genome (NC_000908.2)using minimap2 (v2.24) and samtools (v1.10) to determine the number of reads aligning to the refrence genome:

```
 minimap2 -t8 -ax sr <path/to/reference> <sampleID>_trimmed_R1.fastq.gz <sampleID>_trimmed_R2.fastq.gz | samtools fastq -1 <sampleID>_mapped_R1.fastq.gz -2 <sampleID>_mapped_R2.fastq.gz -f 2
```
Seqeunce reads in both R1 and R2 files for each sample were counted using seqkit stats (v2.3.0):

```
seqkit stats <sampleID>_mapped.R1/R2.fastq.gz> >> <sampleID_outputfile>;done
```
The average read depth was calculated as follows: (number of reads) x (average read length)/length of reference genome

### 3. Percentage of the reference genome covered with 5X read coverage
Reads mapped to the reference genome were aligned to the *M.genitalium* G37 reference genome (NC_000908.2) using snippy(v4.6.0) and samtools (v1.10).  Snippy was also used to determine 

