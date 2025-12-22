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

### 3. Variant calling
Trimmed  paired end reads were aligned to the  *M.genitalium* G37 reference genome (NC_000908.2) using snippy(v4.6.0), requiring a minimum of five or more supporting reads and a variant frequency 0.8 or greater:

```
snippy --mincov 5 --minfrac 0.8 --ref <path/to/reference> --R1 <sampleID>_mapped_R1.fastq.gz  --R2 <sampleID>_mapped_R1.fastq.gz --outdir <sampleID>
```

### 4. Genome coverage
The  percentage of the reference genome covered for each sample was calculated.

Sort snps.bam files found in the snippy output folder for each genome:
```
{#!/bin/bash}
# Set working directory
cd <path to working directory>
# Loop through all subdirectories containing BAM files
for dir in */ ; do
    cd "$dir"
    # Loop through all BAM files in the subdirectory
    for bam_file in *snps.bam; do
        # Get base name of the BAM file
        base=$(basename "$bam_file" .bam)
        # Sort BAM file
        samtools sort -o "$base".sorted.bam "$bam_file"
    done
    cd ..
done
```

To calculate percentage of sites with greater than or equal to 5 times coverage for each snps.sorted.bam file and outputs combined into one csv file:

``` {!/bin/bash}
# Merge the coverage files for all subdirectories and append to the output file
output_file="<path/to/merged_file.csv>"
echo "" > "$output_file"

# Loop through all subdirectories containing BAM files
for dir in */ ; do
    cd "$dir"
    # Loop through all BAM files in the subdirectory
    for bam_file in *.sorted.bam; do
        # Get base name of the BAM file
        base=$(basename "$bam_file" .bam)
        # Calculate genome coverage
        bedtools genomecov -ibam "$bam_file" -d | awk '{if($3>=5) total++}END{if(NR>0) print "'"$dir"' " total/NR*100 "%"}' > "$base"_coverage.txt
        # Append the coverage file for this BAM file to the output file
        cat "$base"_coverage.txt >> "$output_file"
    done
    cd ..
done
```



