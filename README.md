# Separate the bacterial reads using kraken and run the metapro
<https://github.com/ParkinsonLab/MetaPro>

```bash
# Build kraken database
# Building standard database took lot of time and failed at the end.
kraken_stan=https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20251015.tar.gz
output_directory=/fs/scratch/PAS0471/menuka/kraken_stnd_db
ls "$output_directory"
sbatch scripts/build_standard_kraken_db.sh "$kraken_stan" "$output_directory"

## Build the bacterial database fo kraken
bacteria_dir=/fs/scratch/PAS0471/menuka/kraken2_bacteria_db
mkdir $bacteria_dir
cd $bacteria_dir
kraken_cont=oras://community.wave.seqera.io/library/kraken2:2.17.1--1738c34504f3fb18
# this did not work - rsync error
apptainer exec $kraken_cont kraken2-build --download-library bacteria --db bacterial_db
# <https://benlangmead.github.io/aws-indexes/k2> to find diff database
- Add taxonomy file
cp ../kraken_stnd_db/names.dmp .
cp ../kraken_stnd_db/nodes.dmp .
cp ../kraken_stnd_db/seqid2taxid.map .

kraken_bact=https://genome-idx.s3.amazonaws.com/kraken/k2_gtdb_genome_reps_20250609.tar.gz
sbatch /fs/ess/PAS3107/menuka_metatrans/scripts/build_standard_kraken_db.sh "$kraken_bact" "$bacteria_dir"
tail slurm-kraken_database_create-46269504.out

## Run kraken on the bacterial database
outfile=menuka_metatrans/results/kraken_bacteria
mkdir -p "$outfile"/unclassified "$outfile"/classified
sbatch menuka_metatrans/scripts/kraken.sh "$kraken_db" \
    "$outfile" "$report" \
        "$unclassified_output" \
        --cla"$classified_output" \
        --paired "$Read1" "$Read2"
```

# Processing of reads- First focus in viral reads - Virus profiling
**1.** QC: FastQC, MultiQC

```bash
# Run FastQC- to check the quality of sequences
mkdir menuka_metatrans/results/fastqc
outdir=menuka_metatrans/results/fastqc
fastq=$3

for fastq in metaT/data/*.fastq.gz; do
    sbatch menuka_metatrans/scripts/fastqc.sh "$outdir" "$fastq"
done

# Run multiQC to combine the FastQC output
outdir=menuka_metatrans/results/multiqc
sbatch menuka_metatrans/scripts/multiqc.sh \
    menuka_metatrans/results/fastqc \
    menuka_metatrans/results/multiqc
# Lots of sequences were contaminated with adapters.
```
**2.** Trimmomatic and Trimgalore
```bash
## Run Trimmomatic to trim the adapter and bad quality sequences: your samples have lots of primer/adapter sequences
paired=menuka_metatrans/results/trimmomatic/paired
unpaired=menuka_metatrans/results/trimmomatic/unpaired
mkdir $paired $unpaired
for fastq in metaT/data/*1.fastq.gz;do
    R2=${fastq/_1.fastq.gz/_2.fastq.gz} 
    basename_R1=$(basename $fastq) 
    basename_R2=$(basename $R2) 
    sbatch menuka_metatrans/scripts/trimmomatic.sh $fastq $R2 $paired/$basename_R1 \
    $unpaired/$basename_R1 $paired/$basename_R2 $unpaired/$basename_R2
done 

## Look at the reads in which both R1 and R2 survived
grep "Both Surviving" slurm-trim* > trimmomatic_output

## Run FastQC and MultiQC again to see if the quality of read has improved after trimming bad quality sequences
outdir=menuka_metatrans/results/fastqc_trimmomatic
for fastq in menuka_metatrans/results/trimmomatic/paired/*.fastq.gz; do
    sbatch menuka_metatrans/scripts/fastqc.sh "$outdir" "$fastq"
done

sbatch menuka_metatrans/scripts/multiqc.sh \
    menuka_metatrans/results/fastqc_trimmomatic \
    menuka_metatrans/results/multiqc_trimmomatic

# Trimgalore: To remove polyGs, and N
# the sequences still contains lots of Ns and polyGs so remove polyG first using the trimgalore
outdir=menuka_metatrans/results/trimgalore
for R1_fastq in menuka_metatrans/results/trimmomatic/paired/*1.fastq.gz; do
  R2=${R1_fastq/_1.fastq.gz/_2.fastq.gz}
  sbatch menuka_metatrans/scripts/trimgalore.sh "$outdir" "$R1_fastq" "$R2" 
done

# Run MultiQC in trimgalore output
sbatch menuka_metatrans/scripts/multiqc.sh \
    menuka_metatrans/results/trimgalore \
    menuka_metatrans/results/multiqc_trimgalore
```

**3.** Kraken & Krakentools
```bash
# Kraken - to classify and extract the viral reads
ls -lh /fs/scratch/PAS0471/menuka/kraken_stnd_db/
outfile=menuka_metatrans/results/kraken

for R1_fastq in menuka_metatrans/results/trimgalore/*_1.fq.gz; do
  sample_ID=$(basename "$R1_fastq" _1_val_1.fq.gz)
  R2=${R1_fastq/_1_val_1.fq.gz/_2_val_2.fq.gz}
  sbatch menuka_metatrans/scripts/kraken.sh "$outfile/$sample_ID.out" \
  "$outfile/$sample_ID.report" "$R1_fastq" "$R2" 
done

# Extract viral reads using krakentools- 10239 id of viruses
# Clone the krakentool repo
git clone https://github.com/jenniferlu717/KrakenTools.git
KrakenTools/extract_kraken_reads.py --help
# Got an error that No module named 'Biopython'
# Create a conda env and install the dependency Biopython in it.
module load miniconda3/24.1.2-py310
conda create -n krakentools python=3.9 biopython -y
conda activate krakentools
echo $CONDA_PREFIX
# Now this works
KrakenTools/extract_kraken_reads.py --help

## Run KrakenTools
outfile=menuka_metatrans/results/kraken
viral_reads=menuka_metatrans/results/kraken_viral_reads

for R1_fastq in menuka_metatrans/results/trimgalore/*1_val_1.fq.gz; do
  sample_ID=$(basename "$R1_fastq" _1_val_1.fq.gz)
  R2=${R1_fastq/_1_val_1.fq.gz/_2_val_2.fq.gz}
  sbatch menuka_metatrans/scripts/krakentools_viral_reads.sh \
  "$outfile"/"$sample_ID".out \
  "$outfile"/"$sample_ID".report \
  "$R1_fastq" "$R2" \
  "$viral_reads"/"$sample_ID"_1.fastq "$viral_reads"/"$sample_ID"_2.fastq
done 

grep "reads printed to file" slurm* \
> total_viral_reads.out
# Since lots of our viral reads are less than 100k, we will not do the normalization
```

**4.** Metaspades
```bash
# Assemble viral reads to contigs using metaspades, which is part of the spades toolkit
outdir=menuka_metatrans/results/metaspades

for R1 in menuka_metatrans/results/kraken_viral_reads/*_1.fastq; do
    R2=${R1/_1.fastq/_2.fastq}
    sample_ID=$(basename "$R1" _1.fastq)
    sbatch menuka_metatrans/scripts/metaspades.sh "$R1" "$R2" "$outdir"/"$sample_ID"
done

## Copy all the contigs in the same folder and add the sampleID to each of them
for sample_name in menuka_metatrans/results/metaspades/*; do
  sample=$(basename "$sample_name")
  cp "$sample_name"/contigs.fasta \
  menuka_metatrans/results/viral_contigs/"$sample"_contigs.fasta
done

# Use seqkit to filter the contigs less than 300bp
outdir=menuka_metatrans/results/seqkit

for contigs in menuka_metatrans/results/viral_contigs/*.fasta;do
    sample_ID=$(basename "$contigs" )
    sbatch menuka_metatrans/scripts/seqkit.sh "$contigs" "$outdir"/"$sample_ID"
done

# To get the stats of the contigs
apptainer exec $container seqkit stats $outdir/AN_01_contigs.fasta
```

**5.** BlAST
```bash
# Identify which virus the contigs belongs to using BLAST
# Download the refseq of NCBI and blast your contigs against it:https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz, I had viral database in OSC so I will use one
# /fs/ess/PAS0471/condaenv_database_Menuka/NCBI_viral_db/
# Blast our contigs (query) against the viral database(subject)
- BLAST database takes the database prefix, not the directory
db_name=/fs/ess/PAS0471/condaenv_database_Menuka/NCBI_viral_db/refseq_viral
blast=menuka_metatrans/results/blast_viral_contigs

for contigs in menuka_metatrans/results/seqkit/*.fasta; do
sample=$(basename $contigs _contigs.fasta)
sbatch menuka_metatrans/scripts/blast.sh "$contigs" "$db_name" "$blast"/"$sample".tsv
done

# Find out which virus each blast results corresponds to
cut -f2 menuka_metatrans/results/blast_viral_contigs/AN_01.tsv | sort -u >AN_01.txt
head AN_01.txt

# Get the full list of the viral seqeunce present in your database
apptainer exec "$cont" blastdbcmd \
  -db "$db_name" \
  -entry all \
  -outfmt "%a\t%t" > all_viral_names.tsv

# Find out the rota virus presentr in it
grep "rota" all_viral_names.tsv > rota.tsv
# use grep to look for the ID that matches between your sample and database
grep -Ff AN_01.txt all_viral_names.tsv > virus_names.tsv
```

**6.** Bowtie2
```bash
# Run Bowtie2, to map reads to the reference virus genome, this will give counts of reads per virus
# Initial plan was to Use GCA_003087295.1 for RotaC becuase it has more contigsN50 and GCA_002660075.1 for RotaA becuase it was the published in 2017 and have more contigsN50 but later decided to use all Rota genome available in NCBI

A. Activate Conda env to download all the Rata genome from NCBI
datasets download virus genome taxon 28875 --filename rotavirus_data.zip
sbatch menuka_metatrans/scripts/download_genome_NCBI.sh
unzip menuka_metatrans/data/NCBI_rota_genome/rotavirus_data.zip

## Fix the header name of the fasta file its giving an error to create the blast database fo Rota genome
grep ">" ncbi_dataset/data/genomic.fna >ncbi_dataset/data/rota_virus.tsv

sed '/^>/s/ .*//' ncbi_dataset/data/genomic.fna \
> ncbi_dataset/data/new_rotavirus.fna

B. Create a BLAST database
QUERY=ncbi_dataset/data/new_rotavirus.fna
sbatch menuka_metatrans/scripts/blast_db_create.sh "$QUERY"

C. Blast the viral assembly against the NCBI rota database
db_name=menuka_metatrans/results/rota_db/NCBI_rota_genome
out_file=menuka_metatrans/results/blast_NCBI_rota

# Select the results with identity percentage > 80
for fasta in menuka_metatrans/results/viral_contigs/*.fasta;do
sample_ID=$(basename $fasta _contigs.fasta)
sbatch menuka_metatrans/scripts/blast.sh "$fasta" "$db_name" "$out_file/$sample_ID.tsv"
done

# Print only rows with coverage greater than 70% - is it necessary to do it?
for tsv in menuka_metatrans/results/blast_NCBI_rota/*.tsv; do
    sample=$(basename "$tsv")
    awk 'BEGIN{FS=OFS="\t"} $11 >= 70 {print}' "$tsv" \
        > "menuka_metatrans/results/blast_NCBI_rota/$sample.tsv"
done

# Check how many different types of rota virus are present in your blast results
cut -f2 menuka_metatrans/results/blast_NCBI_rota/*.tsv \
  | sort -u > menuka_metatrans/results/blast_NCBI_rota/rotavirus.txt

D. Extract the viral genome from your reference FASTA
# The fasta file is needed to align the reads to the genome
ref_seq=ncbi_dataset/data/genomic.fna # this contains old id
ref_fasta=ncbi_dataset/data/new_rotavirus.fna # this conatins new id
less $ref_fasta

# Use seqkit to extract the fasta of the rota virus accession that was detected in our sample
container=oras://community.wave.seqera.io/library/seqkit:2.13.0--205358a3675c7775
apptainer exec $container seqkit grep -f menuka_metatrans/results/blast_NCBI_rota/rotavirus.txt "$ref_fasta" \
  > menuka_metatrans/results/blast_NCBI_rota/rota_genomes.fna

grep ">" menuka_metatrans/results/blast_NCBI_rota/rota_genomes.fna | wc -l

E. Use Bowtie2 to map the viral reads to the Rota virus genome to get the depth
# Run Bowtie2 to map reads to reference
# First you need to build the Bowtie2 index for the reference genome to efficiently map reads to reference genome
sbatch menuka_metatrans/scripts/bowtie2.sh

# Run Bowtie2 for mapping reads to reference
# Map reads to index - produce BAM : how many reads mapped to each rota virus, then sort/order the alignment by their genomic coordinates position
outdir=menuka_metatrans/results/bowties2_mapping_rota

for R1 in menuka_metatrans/results/kraken_output/kraken_viral_reads/*_1.fastq; do
R2="${R1/_1.fastq/_2.fastq}"
sample_Id=$(basename "$R1" _1.fastq)
sbatch menuka_metatrans/scripts/bowtie2_mapping.sh "$R1" "$R2" "$outdir"/"$sample_Id".sorted.bam
done

F. Count reads per rotavirus 
# First create BAM index file - This allows to quickly find reads aligned in a specific genomic region without scanning whole region
samtools=oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a
for bam in menuka_metatrans/results/bowties2_mapping_rota/*.sorted.bam;do
  apptainer exec "$samtools" samtools index "$bam"
done

mv bowties2.txt menuka_metatrans/results/bowties2_mapping_rota

# Get the counts from the BAM index stats: reference_name,reference_length,mapped_reads, unmapped_reads
for bam in menuka_metatrans/results/bowties2_mapping_rota/*.sorted.bam;do
  sample=$(basename "$bam" .sorted.bam)
  apptainer exec $samtools samtools idxstats "$bam" \
    | awk -v s="$sample" 'BEGIN{OFS="\t"} $3 > 0 {print s,$1,$2,$3}' \
    > menuka_metatrans/results/rotavirus_counts_BAM/$sample.viral_counts.tsv
done

# Combine the counts of each Rota virus mapping to individual sample type in R
# Find out the reads mapping to 
```

```bash
## Run QUAST and BUSCO to pick the best genome out of three rotavirus C genome
# Run quast
outdir=menuka_metatrans/results/quast
for fasta in menuka_metatrans/data/rotavirus_C_reference/GCA*; do
sample_ID=$(basename $fasta _genomic.fna.gz)
sbatch menuka_metatrans/scripts/quast.sh "$fasta" "$outdir"/$sample_ID
done

## Run BUSCO
Get all genomes of Rotavirus
for url in $(cat ftp_links.txt); do
    wget "${url}*_genomic.fna.gz"
done
```
