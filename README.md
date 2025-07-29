# FCS Genome Cleaning Pipeline (with Singularity)

This is a detailed step-by-step pipeline to screen and clean genome assemblies using NCBI's FCS-adaptor and FCS-GX. It includes contaminant masking, overlap analysis with annotations, and downstream functional impact checks.

---

## Step 0: Prepare Environment (e.g., CCAST, ATLAS)

### Log in to your HPC

For CCAST or ATLAS, use:

```bash
ssh your_username@ccast.ndsu.edu      # for CCAST
ssh your_username@atlas-login.ndsu.edu  # for ATLAS
```

> Replace `your_username` with your actual HPC username.

---

### Create a working directory

```bash
mkdir fcs_adapter

cd fcs_adapter
```

### Load required modules

Before running any tools, load the following modules (adjust as needed):

```bash
ml singularity/3.8.0
ml samtools/1.20
ml bedtools/2.30.0
ml blast-plus/2.14.1
ml emboss/6.6.0
```

## 1. Setup & Download Tools

### Load Singularity
```bash
ml singularity/3.8.0
```

### Retrieve FCS-adaptor scripts and image
```bash
curl -LO https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh
chmod 755 run_fcsadaptor.sh

curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-adaptor.sif -Lo fcs-adaptor.sif
```

---

## 2. Screen Genome with `run_fcsadaptor.sh`

```bash
mkdir fcs_adapter_cleaned

./run_fcsadaptor.sh --fasta-input /directory/where/you/saved/your/wheat.fasta \
  --output-dir ./fcs_adapter_cleaned \
  --euk \
  --container-engine singularity \
  --image fcs-adaptor.sif
```

### Example PBS Script for Cluster
```bash
#!/bin/bash
#PBS -q bigmem
#PBS -N FCS_adaptor_screen_cleaned
#PBS -l select=1:mem=512gb:ncpus=48
#PBS -l walltime=72:00:00
#PBS -M user.name@ndsu.edu
#PBS -m abe
#PBS -W group_list=x-ccast-prj-name

ml singularity/3.8.0
cd /directory/where/this/saved/fcs_adapter

./run_fcsadaptor.sh --fasta-input /directory/where/you/saved/your/wheat.fasta \
  --output-dir /directory/where/this/saved/fcs_adapter_cleaned \
  --euk --container-engine singularity --image fcs-adaptor.sif
```

---

## 3. Clean Genome Using FCS-GX

### Download scripts and container
```bash
curl -LO https://github.com/ncbi/fcs/raw/main/dist/fcs.py
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-gx.sif -Lo fcs-gx.sif
export FCS_DEFAULT_IMAGE=fcs-gx.sif
```

### Run Cleaning
```bash
ml singularity/3.8.0

python3 ./fcs.py clean genome \
  --input /directory/where/this/saved/wheat.fasta \
  --action-report/directory/where/this/saved/fcs_adapter_cleaned/fcs_adaptor_report.txt \
  --output /directory/where/this/saved/fcs_adapter_cleaned/wheat_masked.fasta \
  --contam-fasta-out /directory/where/this/saved/fcs_adapter_cleaned/contam_sequences.fasta
```

---

## 4. Mask Genome Instead of Trimming

```bash
sed -i 's/ACTION_TRIM/FIX/g' fcs_adaptor_report.txt
```

### Re-run FCS-GX with FIX applied
```bash
python3 ./fcs.py clean genome \
  --input /directory/where/this/saved/wheat.fasta \
  --action-report /directory/where/this/saved/fcs_adapter_cleaned_wheat/fcs_adaptor_report.txt \
  --output /directory/where/this/saved/fcs_adapter_cleaned_wheat/wheat_masked.fasta \
  --contam-fasta-out /directory/where/this/saved/fcs_adapter_cleaned/contam_sequences.fasta
```

---

## 5. Extract and Summarize Contaminated Regions

### Extract lengths
```bash
awk '
BEGIN { FS="\t"; OFS="\t" }
{
    for (i=1; i<=NF; i++) {
        if ($i ~ /[0-9]+\.\.[0-9]+/) {
            range_col = i;
            break;
        }
    }
    if (range_col) {
        split($range_col, ranges, ",");
        for (i in ranges) {
            if (match(ranges[i], /([0-9]+)\.\.([0-9]+)/, arr)) {
                print $1, arr[1], arr[2], arr[2] - arr[1] + 1;
            }
        }
    }
}' fcs_adaptor_report.txt > masked_wheat_lengths
```

### Get totals
```bash
awk '{sum += $4} END {print "Total:", sum}' masked_wheat_lengths
awk '/^>/ {next} {total += length($0)} END {print "Total bases:", total}' wheat_masked.fasta
```

---

## 6. Generate BED for Contaminants and Extract Sequences

### Create BED file
```bash
awk '
/FIX/ {
    gsub(/>/, "", $0);
    split($4, ranges, ",");
    for (i in ranges) {
        if (match(ranges[i], /([0-9]+)\.\.([0-9]+)/, coords)) {
            start = coords[1] - 1;
            end = coords[2];
            print $1, start, end;
        }
    }
}' fcs_adaptor_report.txt > contaminated_regions.bed

awk '{$1=$1; OFS="\t"; print}' contaminated_regions.bed > wheat_contaminated_regions_tab.bed
```

### Index and extract
```bash
ml samtools/1.20
samtools faidx wheat.fasta

ml bedtools/2.30.0
bedtools getfasta -fi wheat.fasta -bed wheat_contaminated_regions_tab.bed -fo wheat_contaminated_sequences.fasta
```

---

## 7. Intersect Contaminants with Annotations

### Generate Annotation Feature Files
```bash
awk '$3 == "exon"' wheat.gff3 > wheat_exons.gff3
awk '$3 == "CDS"' wheat.gff3 > wheat_cds.gff3
awk '$3 == "five_prime_UTR" || $3 == "three_prime_UTR"' wheat.gff3 > wheat_utrs.gff3
awk '$3 == "gene"' wheat.gff3 > wheat_genes.gff3
```

### Get introns from BED subtraction
```bash
awk '{OFS="\t"; print $1, $4-1, $5, $9, ".", $7}' wheat_exons.gff3 > wheat_exons.bed
awk '{OFS="\t"; print $1, $4-1, $5, $9, ".", $7}' wheat_genes.gff3 > wheat_genes.bed

bedtools subtract -a wheat_genes.bed -b wheat_exons.bed > wheat_introns.bed
```

### Intersections
```bash
bedtools intersect -a wheat_contaminated_regions_tab.bed -b wheat_exons.bed -wa -wb > wheat_contaminated_in_exons.txt
bedtools intersect -a wheat_contaminated_regions_tab.bed -b wheat_introns.bed -wa -wb > wheat_contaminated_in_introns.txt
bedtools intersect -a wheat_contaminated_regions_tab.bed -b wheat_cds.gff3 -wa -wb > wheat_contaminated_in_cds.txt
bedtools intersect -a wheat_contaminated_regions_tab.bed -b wheat_utrs.gff3 -wa -wb > wheat_contaminated_in_utrs.txt
bedtools intersect -a wheat_contaminated_regions_tab.bed -b wheat_genes.bed -v > wheat_contaminated_intergenic.txt
```

### Summary
```bash
echo "Exons: $(wc -l < wheat_contaminated_in_exons.txt)"
echo "Introns: $(wc -l < wheat_contaminated_in_introns.txt)"
echo "CDS: $(wc -l < wheat_contaminated_in_cds.txt)"
echo "UTRs: $(wc -l < wheat_contaminated_in_utrs.txt)"
echo "Intergenic: $(wc -l < wheat_contaminated_intergenic.txt)"
```

---

## 8. Extract and Compare Protein Sequences

```bash
awk '/^>TRAES.wheat.r1.4BG00763140.1/{print; found=1; next} /^>/ {found=0} found' TRAES.wheat.pgsb.r1.Mar2024.aa.fa > wheat_extracted_sequence.fasta
```

---

## 9. Translation and BLAST

### Translate sequences
```bash
ml emboss/6.6.0
transeq -sequence cds85.fasta -outseq protein_original.fasta
transeq -sequence contaminated_sequence.fasta -outseq protein_contaminated.fasta
```

### BLAST alignment
```bash
ml blast-plus/2.14.1
makeblastdb -in protein_original.fasta -dbtype prot -out protein_db
blastp -query protein_contaminated.fasta -db protein_db -out alignment_results.txt -outfmt 6
```

---

## 10. Additional Gene/Contaminant Extraction

```bash
awk -F 'ID=' '/gene/ {print $2}' wheat_affected_annotations.gff3 | cut -d';' -f1 | sort | uniq > wheat_gene_ids.txt

awk '$3 == "gene" {print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' wheat_affected_annotations.gff3 > wheat_genes.bed
bedtools getfasta -fi wheat_pm_v1.fasta -bed wheat_genes.bed -fo wheat_sequences.fasta

makeblastdb -in wheat_sequences.fasta -dbtype nucl -out wheat_db

ml blast/2.16.0
blastn -query wheat_gene_ids.txt -db nt -out wheat_blast_results.txt -evalue 1e-5 -outfmt 6
```

---

## Maintainer

This guide was written by **Ruby Mijan** as a user-friendly walkthrough for applying NCBI’s FCS tools.

- Tools developed by: [NCBI FCS Team](https://github.com/ncbi/fcs)
> This documentation only organizes and demonstrates how to apply the tools using real genome data in an HPC workflow.
