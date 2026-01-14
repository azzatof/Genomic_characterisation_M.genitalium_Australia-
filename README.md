# Genomic characterisation of *Mycoplasma genitalium* in Victoria, Australia, reveals lineage diversification and drivers of antimicrobial resistance
Francesca Azzato, George Taiaroa , Janath Fernando, Mona L Taouk, Vesna De Petra, Lenka A Vodstrcil, Erica L Plummer, Kerry Raios, Danielle J Ingle, Niamh Meagher, Jacqueline Prestedge, Eloise Williams, Leon Caly, Benjamin P Howden, Shivani Pasricha, Catriona S Bradshaw, Deborah A Williamson

This GitHub repository contains all code used in this study. Sequence reads are available from the NCBI database under BioProjects PRJNA1367946 and PRJEB5172, with accession numbers provided in the Supplementary Dataset accompanying this publication. All analyses performed in this study can be replicated using the code provided, with input names, directory paths, and file names modified as required.


## Quality Control of *Mycoplasma genitalium* sequence reads

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
### 4. Kraken
Mapped paired end reads for each genome were input into kraken2 (v2.12 using the Plus PF database):

```
kraken2 --db /path/to/KRAKEN2/PLUSPF --gzip-compressed --paired <path/to/smapleID_mapped_R1.fastq.gz> <path/to/sampleID_mapped_R2.fastq.gz> --report <sampleID_report_.txt>
```
Top species match and percentage were extrated from each result file and summarised:
```
 for file in *.txt; do awk 'FNR==1 || $1>=16 {print FILENAME, $0}' "$file"; done > merged_Kraken_topID.csv
```

## Assemblies
### 1. Shovill
*De novo* genome assemblies were generated using Shovill (v1.1.0)

```
shovill --gsize 0.58M -- outdir <sampleID> -- R1 <sampleID>_trimmed_R1.fastq.gz  --R2 <sampleID>_trimmed_R2.fastq.gz
```
### 2. Quality control of genome assemblies

To assess the quality of  *De novo* genome assemblies QUAST (v5.3.0) was performed.

```
{#!/bin/bash}
# Set working directory
cd path/to/assemblies/contigs.fa
#loop through all subdirectories
for dir in */; do
    cd "$dir" || continue
    # Loop through all contig.fa files in subdirectory
    for file in contigs.fa; do
        if [[-f "$file" ]]; then
           echo "Running Quast on $dir$file"
           quast "$file"
        fi
    done
```
 
## Whole genome Maximum likelihood phylogeny (ML)

### 1. Alignment of Genomes to *M.genitalium* reference genome
A pseudoalignment of all genomes was generated using Snippy (v4.6.0), based on the variant calling output described in the quality control section above.

The snippy-core function in Snippy (v4.6.0) was used to align all genomes against the Mycoplasma genitalium G37 reference genome (GenBank accession number NC_000908.2):

```
snippy-core --ref G37.fasta Snippy/output/folder/for/each/genome
```

### 2. Recombination filtering
High recombination regions of the MgpA operon as identifed previously by Fookes, *et al*. were masked to N characters in the snippy full alignment using BEDTools (v2.3.0).  These problematic positions are listed in a bed file.
```
bedtools maskfasta -fi MG.core.full.aln -fo MG.core.full.positional.aln -bed Global_local_MgpA_reco_sites.bed
```

Additional recombinant regions to the above sites were identified using Gubbins (v2.4.1), appying a minumum SNP threshold of 40 and using the maksed alignment generated in the above step:

```
run_gubbins.py -i 20 -- threads 10 min_snps 40 MG.core.full.positional.aln
```
Following Gubbins analysis, a SNP alignment comprising 7,097 variable sites was generated using Core-SNP-filter (v0.1.1). A 95% soft-core threshold was applied, using the Gubbins polymorphic sites alignment as input.

```
coresnpfilter -c 0.95 MG.core.full.positional.aln > MG.Alignment.positional.filtered_polymorphic_sites_95.fasta
```

### 3. Maximum likelihood phylogeny
A maximum likelihood phylogenetic tree was generated using IQ-tree (v2.2.0.3):

```
iqtree -s MG.Alignment.positional.filtered_polymorphic_sites_95.fasta -B 1000
```

## Population Structure 
### 1. Heirarchial Bayesian of Population Sturcture (BAPS)
Lineages were identified using rhierbaps (v1.1.4) algorithm in R using a maximum depth of three. Up to 50 populations were considered and assignment probabilities were taken into considertion.

The following R code was used:
```
library("ggtree")
library("phytools")
library("rhierbaps")

set.seed (1234)

MG.core.full.aln <- read.dna('MG.Alignment.positional.filtered_polymorphic_sites_95.fasta', format="fasta" )
snp.matrix <- load_fasta("MG.Alignment.positional.filtered_polymorphic_sites_95.fasta")

hb.results <- hierBAPS(snp.matrix, max.depth = 3, n.pops = 50, quiet = TRUE)
head(hb.results$partition.df)

write.csv(hb.results$partition.df, file = "hierbaps_partition.csv", col.names = TRUE, row.names = FALSE)
```
Results:

- [BAPS result table](results/hierbaps_partition.csv)

## Ancestral reconstruction
To assess the evolutionary trajectory of *M.genitalium* ancestral character state reconstruction was performed using ape (v4.3.3) in R under the equal rates model.  The maximum likelihood tree was used as input.

The following R code was used
```
library(ape)  # use for reading tree
library(phytools)  # use for ASR
library(tidyverse)  # data manipulation
library(ggplot2)  # plotting
library(ggtree)  # tree plotting


# Read in the tree file
tree <- read.tree("Alignment.positional.filtered_polymorphic_sites_95.fasta_midpoint.nwk")

# Load metadata
traits <- read.csv("meta_data_for_rectangular_tree.csv", strip.white = TRUE, header = TRUE) %>%
  select(Isolate, BAPS) %>%
  rename(tips = Isolate, lineage = BAPS)

# Set names as tip labels
trait_vector <- setNames(traits$lineage, traits$tip)

# This should return TRUE if all tip names in the vector match those in the tree
all(tree$tip.label %in% names(trait_vector))

# Run ancestral state (default model= ER)
asr_result <- ace(trait_vector, tree, type = "discrete")

# View likelihood of each BAP being that "state" for each node
head(asr_result$lik.anc)

# Make that a dataframe
asr_df <- as.data.frame(asr_result$lik.anc)
asr_df$node <- as.numeric(rownames(asr_df))
rownames(asr_df) <- NULL

# Reshape to long format
df_long <- asr_df %>%
  pivot_longer(cols = starts_with("BAPS"),
               names_to = "BAPS_group",
               values_to = "likelihood")

# For each node, select the BAPS group with the highest likelihood
df_max <- df_long %>%
  group_by(node) %>%
  slice_max(likelihood, with_ties = FALSE) %>%
  ungroup()
```
## cgMLST
### 1. Create Schema
The *M.genitalium* cgMLST schema was developed in this study.
A prodigal training file for *M.genitalium* was made using the *M.genitalium* reference genome ((G37; Genbank accession NC_000908.2), using prodigal (v2.6.3):

```
prodigal -i G37.fna -t MG_training_file.trn -p single
```
Training file: [MG_training_file](files/MG_training_file/trn)

The cgMLST scheme was developed using the Create Schema module in chewBBACA (version 3.3.2):

```
chewBBACA.py CreateSchema -i path/to/input/assemblies/folder -g path/to/output/folder --n MG_schema --ptf path/to/MG/Training/file/MG_training_file.trn
```

### 2. Annotate Schema
Gene annotation for the loci included in the schema was performed using the UniProtFinder module:

```
chewBBACA.py Uniprotfinder -i path/to/schmea/folder/MG_cgmlst_Schema -o path/to/output/folder/schema_annotations -t /path/to/cds_coordinates.tsv --taxa "Mycoplasma genitalium" --cpu 4
```
### 3. Allele Calling
Allele calling for all *M. genitalium* genomes was peformed using the AlleleCall module, with the genome assemblies provided as input:

```
chewBBACA.py AlleleCall -i path/to/assemblies/input/folder -g path/to/schmea/folder/MG_cgmlst_Schema -o path/to/output/folder/MG_cgmlst_allele_call_2 --cpu 4
```
### 4.Define the set of loci constituting the dataset’s core genome (ExtractCgMLST)
To determine the loci that are  present in more than 95% of the dataset ExtractCgMLST module was used:

```
chewBBACA.py ExtractCgMLST -i path/to/input/file/MG_cgmlst_allele_call_2/results_alleles.tsv  -o path/to/output/folder/MG_extracted_cgmlst_2
```
cgMLST allele calling final results:
Table of alleles for each gene in the schema for all isolates: [cgMLST95.tsv](results/cgMLST/cgMLST95.tsv)

1760 loci included in the final schema for this study: [cgMLSTschema95.txt](results/cgMLST/cgMLSTschema95.txt)

Number of and percentage of uncalled genes per genome in this dataset: [mdata_stats.tsv](results/cgMLST/mdata_stats.tsv)

Gene presence/absence table (1 = present. 0 = absent): [Presence_Absence.tsv](results/cgMLST/Presence_Absence.tsv)

### 5. cgMLST allelic differences
The cgMLST results table was converted to a symmetrical distance matrix using cgMLST-dists:

```
cgMLST-dists cgMLST95.tsv > cgmlst_distances_0.95.txt
```
Result matrix: [cgmlst_distances_0.95.txt](results/cgMLST/cgmlst_distances_0.95.txt)


## Adjusted Rand Index
To assess the stability of BAPS group assignments we calculated the adjusted rand index (ARI) between the original BAPS group clustering assignmnets  and those derived from 1000 botstrap and subsampled data sets of 50%, 70% and 90%.  To calculate ARI we used mclust (v6.1.1).

### 1. ARI across bootstrap replicates

To calculate the distribution of ARI scores across 1,000 bootstrap iterations of the hierarchical Bayesian Analysis of Population Structure (BAPS), the following R code was used:

```
library(phangorn)
library("mclust")  # for adjustedRandIndex
library(rhierbaps)



# Set seed for reproducibility
set.seed(1234)

# Load SNP alignment
aln_file <- "Alignment.positional.filtered_polymorphic_sites_95.fasta"
snp.matrix <- load_fasta(aln_file)

# Run original hierBAPS
original_result <- hierBAPS(snp.matrix, max.depth = 3, n.pops = 20, quiet = TRUE)

# Extract original clustering (level 1)
original_partition <- original_result$partition.df
orig_clusters <- original_partition$cluster1

# Function to perform one bootstrap replicate
bootstrap_hierbaps <- function(snp_matrix, original_clusters) {
  n_cols <- ncol(snp_matrix)
  
  # Resample SNP positions with replacement
  resampled_snp <- snp_matrix[, sample(1:n_cols, replace = TRUE)]
  
  # Run hierBAPS
  boot_result <- hierBAPS(resampled_snp, max.depth = 3, n.pops = 20, quiet = TRUE)
  
  # Extract level 1 clustering
  boot_partition <- boot_result$partition.df
  
  # Ensure sample order matches
  boot_clusters <- boot_partition$cluster1[match(rownames(original_partition), rownames(boot_partition))]
  
  # Compute Adjusted Rand Index
  ari <- adjustedRandIndex(original_clusters, boot_clusters)
  return(ari)
}

# Run bootstraps
n_bootstraps <- 1000
bootstrap_ARI <- numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  cat("Running bootstrap", i, "\n")
  bootstrap_ARI[i] <- bootstrap_hierbaps(snp.matrix, orig_clusters)
}

# Results
summary(bootstrap_ARI)
hist(bootstrap_ARI, breaks = 20, main = "Bootstrap Adjusted Rand Index", xlab = "ARI")

# Save results
write.csv(bootstrap_ARI, "bootstrap_ARI_results.csv", row.names = FALSE)
```

Results: 
Distribution of ARI scores for 1000 bootstrap replicates: [bootstrap_ari_scores_1000.csv](results/Stats/bootstrap_ari_scores_1000.csv)


## 2. ARI across subsampling fractions
To calculate the distribution of ARI scores across subsampling fractions for hierarchical Bayesian Analysis of Population Structure (BAPS), the following R code was used:

```
# Load packages
library(rhierbaps)
library(ape)
library(dplyr)
library(mclust)
library(ggplot2)
library("Hmisc")
setwd("~/Desktop/tree_analysis/MG_AMR_graphs/SNP distance threshold")


# ----- Parameters -----
alignment_path <- "Alignment.positional.filtered_polymorphic_sites_95.fasta"  # multi-FASTA alignment used in rhierbaps
metadata_path <- "fullclusters.csv"           # must include genome IDs and original rhierbaps cluster
subsample_fractions <- c(0.9, 0.7, 0.5)  # Subsampling fractions to test: 90%, 70%, 50%
n_replicates <- 10  # Number of replicates for each subsampling fraction
ari_results <- data.frame(SubsampleFraction = numeric(0), ARI = numeric(0))  # To store ARI results
output_dir <- "rbaps_robustness_output"  # Output directory
dir.create(output_dir, showWarnings = FALSE)  # Create the output directory

# ----- Load Data -----
alignment <- load_fasta(alignment_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Make sure genome names in metadata match alignment
metadata <- metadata %>% filter(Isolate %in% rownames(alignment))
alignment <- alignment[metadata$Isolate, ]  # re-order

# Store original cluster labels
full_labels <- metadata$Cluster  # should be a factor or integer

# ----- Stratified Subsampling Function -----
stratified_subsample <- function(metadata, group_col, subsample_fraction) {
  metadata %>%
    group_by(!!sym(group_col)) %>%
    sample_frac(size = subsample_fraction) %>%
    ungroup()
}

# ----- Function to Run rhierbaps and Calculate ARI -----
run_baps_subsample <- function(sub_metadata, full_labels, alignment) {
  # Subset alignment
  sub_genomes <- sub_metadata$Isolate
  sub_alignment <- alignment[sub_genomes, ]
  
  # Run rhierbaps (1 level of clustering)
  baps_result <- hierBAPS(sub_alignment, max.depth = 2, n.pops = 20)
  
  # Extract cluster labels (level 1)
  baps_clusters <- baps_result$partition.df
  sub_labels <- baps_clusters %>%
    filter(Isolate %in% sub_genomes) %>%
    arrange(match(Isolate, sub_genomes)) %>%
    pull(`level 1`)
  
  # Match full_labels for these genomes
  ref_labels <- full_labels[match(sub_genomes, metadata$Isolate)]
  
  # Check if lengths are the same
  if (length(ref_labels) != length(sub_labels)) {
    stop(paste("Length mismatch: ref_labels =", length(ref_labels), 
               "sub_labels =", length(sub_labels)))
  }
  
  # Compute ARI
  ari <- adjustedRandIndex(ref_labels, sub_labels)
  
  return(ari)
}

# ----- Run Subsampling and BAPS for Multiple Fractions -----
all_ari_scores <- data.frame(SubsampleFraction = numeric(0), replicate = numeric(0), ARI = numeric(0))

# Loop over each subsampling fraction (90%, 70%, 50%)
for (subsample_fraction in subsample_fractions) {
  cat("Running subsample fraction:", subsample_fraction, "\n")
  
  # Store ARI scores for this subsampling fraction
  ari_scores <- numeric(n_replicates)
  
  # Run replicates
  set.seed(123)  # reproducibility
  for (i in 1:n_replicates) {
    sub_metadata <- stratified_subsample(metadata, "Cluster", subsample_fraction)
    ari <- tryCatch({
      run_baps_subsample(sub_metadata, full_labels, alignment)
    }, error = function(e) {
      warning(paste("Replicate", i, "failed:", e$message))
      return(NA)
    })
    
    ari_scores[i] <- ari
  }
  
  # Store results for this subsample fraction
  all_ari_scores <- rbind(all_ari_scores, data.frame(SubsampleFraction = subsample_fraction, replicate = 1:n_replicates, ARI = ari_scores))
}

# ----- Save and Plot Results -----
# Save results to CSV
write.csv(all_ari_scores, file.path(output_dir, "ari_scores.csv"), row.names = FALSE)
```
Results: 
Distribution of ARI scores for subsampling fractions: [ARI_subsampling_fractions_results.csv](results/Stats/ARI_subsampling_fractions_results.csv)

### Linear Regression 
To assess the relationship between pairwise SNP distances and pairwise allelic distances we applied linear regression model.

The following r-code was used with [

```
dist <- read.csv("cgmlst_dist_melt_final.csv",header=TRUE)
