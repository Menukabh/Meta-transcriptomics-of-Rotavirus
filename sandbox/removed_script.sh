# Run Krakentools to extract just the viral reads - Run for one sample first
outfile=menuka_metatrans/results/kraken/AN_01.out
report=menuka_metatrans/results/kraken/AN_01.report
Read1=menuka_metatrans/results/trimgalore/AN_01_1_val_1.fq.gz
Read2=menuka_metatrans/results/trimgalore/AN_01_2_val_2.fq.gz
viral_reads1=menuka_metatrans/results/kraken_viral_reads/AN_01_1_val_1.fq.gz
viral_reads2=menuka_metatrans/results/kraken_viral_reads/AN_01_2_val_2.fq.gz
sbatch menuka_metatrans/scripts/viral_reads.sh "$outfile" \
  "$report" \
  "$Read1" \
  "$Read2" \
  "$viral_reads1" \
  "$viral_reads2"

# Check whether the NCBI reference database contain the virus infromation and when was it built
module spider blast
module load blast-database/2024-07
cont=oras://community.wave.seqera.io/library/blast:2.17.0--3d1eb1104ccfd59c

# To get some info- dtabase name, total sequences, bases
apptainer exec $cont blastdbcmd \
  -db /fs/ess/PAS0471/condaenv_database_Menuka/NCBI_viral_db/refseq_viral \
  -info
  
# Check the header
apptainer exec $cont blastdbcmd \
  -db /fs/ess/PAS0471/condaenv_database_Menuka/NCBI_viral_db/refseq_viral \
  -entry all | head

ls -lh --time-style=long-iso \
  /fs/ess/PAS0471/condaenv_database_Menuka/NCBI_viral_db/refseq_viral.nsq


# I had .gz extension in my file but it was not zipped so remove it from all files name
for fastq in menuka_metatrans/results/kraken_viral_reads/*.fastq.gz; do
  mv "$fastq" "${fastq%.gz}"
done

# Download porcine rotavirusC to check the quality of the genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/087/275/GCA_003087275.1_ASM308727v1/GCA_003087275.1_ASM308727v1_genomic.fna.gz \
-P menuka_metatrans/data/rotavirus_C_reference


# Manipulating Blast results of the viral contigs
cut -f2 menuka_metatrans/results/blast_viral_contigs/AN_01.tsv \
| sort -u | wc -l 
cut -f2 menuka_metatrans/results/blast_NCBI_rota/* \
  | sort -u > menuka_metatrans/results/blast_NCBI_rota/rotavirus.txt

wc -l menuka_metatrans/results/blast_NCBI_rota/rotavirus.txt

## Select the blast results with the identifty greater than 80%
awk 'BEGIN{FS=OFS="\t"} $3 >= 80 {print $2}' \
  menuka_metatrans/results/blast_viral_contigs/*.tsv \
  | sort -u > menuka_metatrans/blast/viral_accessions.txt
wc -l menuka_metatrans/blast/viral_accessions.txt

# Extract the viral genome from your reference FASTA
ref_fasta=/fs/ess/PAS0471/condaenv_database_Menuka/NCBI_viral_db/refseq_viral.fna
less $ref_fasta

# Use seqkit to extract the fasta for the accession that
container=oras://community.wave.seqera.io/library/seqkit:2.13.0--205358a3675c7775
apptainer exec $container seqkit grep -f menuka_metatrans/blast/viral_accessions.txt "$ref_fasta" \
  > menuka_metatrans/blast/selected_viral_genomes.fna

# Build the index for the bowtie2
sbatch menuka_metatrans/scripts/bowtie2.sh

## Map reads to index - produce BAM : how many reads mapped to each MAGs, 
# then sort/order the alignment by their genomic coordinates position
outdir=menuka_metatrans/results/bowties2_mapping
for R1 in menuka_metatrans/results/kraken_output/kraken_viral_reads/*.fastq; do
R2="${R1/_1.fastq/_2.fastq}"
sample_Id=$(basename "$R1" _1.fastq)
sbatch menuka_metatrans/scripts/bowtie2_mapping.sh "$R1" "$R2" "$outdir"/"$sample_Id".sorted.bam
done

## Count reads per virus/reference
# First create the index file
samtools=oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a
for bam in menuka_metatrans/results/bowties2_mapping/*.sorted.bam;do
  apptainer exec "$samtools" samtools index "$bam"
done

# Get the counts from the BAM index stats: reference_name,reference_length,mapped_reads, unmapped_reads
for bam in menuka_metatrans/results/bowties2_mapping/*.sorted.bam;do
  sample=$(basename "$bam" .sorted.bam)
  apptainer exec $samtools samtools idxstats "$bam" \
    | awk -v s="$sample" 'BEGIN{OFS="\t"} $3 > 0 {print s,$1,$2,$3}' \
    > menuka_metatrans/results/viral_counts_BAM/${sample}.viral_counts.tsv
done

## combine all samples
cat menuka_metatrans/results/viral_counts_BAM/*viral_counts.tsv \
  > all_samples_viral_counts.tsv

mkdir -p menuka_metatrans/results/merged_virome

for blast in menuka_metatrans/results/blast_viral_contigs/*.tsv
do
  sample=$(basename "$blast" .tsv)

  awk -v s="$sample" 'BEGIN{FS=OFS="\t"}
  {
    virus=$2
    contigs[virus][$1]=1
    if ($12 > best_bitscore[virus]) best_bitscore[virus]=$12
    if ($3 > best_pident[virus]) best_pident[virus]=$3
    if ($13 > best_qcov[virus]) best_qcov[virus]=$13
  }
  END{
    for (v in contigs) {
      n=0
      for (c in contigs[v]) n++
      print s,v,n,best_pident[v],best_qcov[v],best_bitscore[v]
    }
  }' "$blast"

done > menuka_metatrans/results/merged_virome/blast_summary_by_sample_virus.tsv

# Merge blast summary with read mapping counts
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
  key=$1 FS $2
  genome_len[key]=$3
  mapped_reads[key]=$4
  next
}
{
  key=$1 FS $2
  print $0, genome_len[key], mapped_reads[key]+0
}' \
all_samples_viral_counts.tsv \
menuka_metatrans/results/merged_virome/blast_summary_by_sample_virus.tsv \
> menuka_metatrans/results/merged_virome/final_virome_profile.tsv


# I was getting some issue with the blastdatabase, so checked the integrity of the database
# Check the database
module load blast-plus/2.16.0
blastdbcheck -db $db_name -dbtype nucl
blastdbcmd -db $db_name -info
# Header of database
blastdbcmd -db $db_name -entry all | head
blastdbcmd -db $db_name -entry all -outfmt "%i" | head
blastdbcmd -db $db_name -entry all | head -50
blastdbcmd -db $db_name -entry all | grep "-"
grep -v "^>" query.fasta | grep -n "[^ACGTNacgtn]"



