# Merge count data generated from BAM file with the blast data
# Load packages:
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(writexl)
library(tidyverse)

# Files to combine
virus_file  <- "ncbi_dataset/data/rota_virus.tsv" # This file has all the description
count_dir <- "menuka_metatrans/results/rotavirus_counts_BAM"

# Read and combine all count files
count_files <- list.files(count_dir,
  pattern = "viral_counts.tsv",
  full.names = TRUE)

# Read count data with four columns
counts <- purrr::map_dfr(
  count_files,
  readr::read_tsv,
  col_names = c("sample", "accession", "reference_length", "mapped_reads"),
  col_types = "ccnn")

# read virus info when accession and other description in other columns
lines <- readLines(virus_file)
head(lines)

# Separate accession name with the description using sub command
virus_info <- data.frame(
  accession = sub("^>(\\S+).*", "\\1", lines),
  description = sub("^>\\S+\\s+", "", lines),
  stringsAsFactors = FALSE)

virus_info$rotavirus <- sub(".*(([A-Za-z]+ )?[Rr]otavirus( [A-Z])?).*", 
                            "\\1", virus_info$description)
virus_info$isolate <- sub(".*(isolate [^ ]+).*", "\\1", virus_info$description)
virus_info$genotype <- sub(".*(VP[0-9]+/G[0-9A-Z]+|G[0-9A-Z]+|P\\[[0-9]+\\]).*", "\\1", virus_info$description)

virus_info$gene <- NA
virus_info$gene[grepl("VP7", virus_info$description)] <- "VP7"
virus_info$gene[grepl("VP4", virus_info$description)] <- "VP4"
virus_info$gene[grepl("VP6", virus_info$description)] <- "VP6"
virus_info$gene[grepl("NSP1", virus_info$description)] <- "NSP1"
virus_info$gene[grepl("NSP2", virus_info$description)] <- "NSP2"
virus_info$gene[grepl("NSP3", virus_info$description)] <- "NSP3"
virus_info$gene[grepl("NSP4", virus_info$description)] <- "NSP4"
virus_info$gene[grepl("NSP5", virus_info$description)] <- "NSP5"


#virus_info$accession <- gsub("^>", "", virus_info$accession)
merged <- counts %>%
  left_join(virus_info, by = "accession")

write_xlsx(merged, "rota_read_mapping_refseq.xlsx")

# Select the ones that belongs to porcine
porcine <- merged %>%
  filter(grepl("Porcine", description))

#################################################################
# Separate the multiple genotypes G, P's and I
data <- read_excel("rota_count_blast.xlsx")

data$G <- sub("([a-zA-Z]+)([0-9]+).*", 
              "\\1\\2", data$genotype)
data$P <- str_extract(data$genotype, "P[0-9]+")
data$I <- str_extract(data$genotype, "I[0-9]+")
data <- data |> 
  mutate(gvalue = str_extract_all(data$genotype, "G\\d+"))
data$G2 <- sapply(b$gvalue, `[`, 2)
data$G3 <- sapply(b$gvalue, `[`, 3)
data <- data |> 
  select(!gvalue)

write_xlsx(data, "rota_genotype_separated.xlsx")

#################################################################
# Work with the blast data
blast_dir <- "menuka_metatrans/results/blast_NCBI_rota"
blast_files <- list.files(blast_dir,
                          pattern = ".tsv",
                          full.names = TRUE)

blast_data <- purrr::map_dfr(
  blast_files,
  ~readr::read_tsv(.x,
  col_names = c("qseqid", "accession", "pident", "length", "mismatch", "gapopen", "qstart", "qend",  "sstart", "send", "evalue", "bitscore", "qcovs"),
  col_types = "ccnnnnnnnnnnn") |> 
  mutate(sample_id = tools::file_path_sans_ext(basename(.x))))
length(unique(blast_data$sample_id))

# Filter data with the alignment length > 700bp
blast_700bp <- blast_data |> 
  filter(length > 700) 

blast_700_cov <- blast_700bp |> 
  filter(qcovs > 80) 

length(unique(blast_700bp$sample_id))
length(unique(blast_700_cov$sample_id))

# Combine this data with the virus name
blast_hits <- blast_700_cov %>%
  left_join(virus_info, by = "accession")

unique(blast_hits$rotavirus)

# Select only the pig or porcine rota hits
porcine_rota <- blast_hits |> 
  filter(grepl("porcine|pig", description, ignore.case = TRUE))

length(unique(porcine_rota$sample_id))

write_xlsx(porcine_rota, "blast_hits_summary.xlsx")

##########################################################
#################################################################
# Separate the multiple genotypes G, P's and I
data <- read_excel("blast_hits_summary.xlsx")

data$G <- sub("([a-zA-Z]+)([0-9]+).*", 
              "\\1\\2", data$genotype)
data$P <- str_extract(data$genotype, "P[0-9]+")
data$I <- str_extract(data$genotype, "I[0-9]+")
write_xlsx(data, "rota_blast_genotype_separated.xlsx")

data <- data |> 
  mutate(gvalue = str_extract_all(data$genotype, "G\\d+"))
data$G2 <- sapply(data$gvalue, `[`, 2)

data$G2 <- sapply(b$gvalue, `[`, 2)
data$G3 <- sapply(b$gvalue, `[`, 3)
data <- data |> 
  select(!gvalue)

write_xlsx(data, "rota_blast_genotype_separated.xlsx")
  
###################################################
# Now group by samples and see what are the different types of rota each sample has
blast_hits_rota_count <- data |> 
  group_by(sample_id) |>  
  summarise(n_unique_rota = n_distinct(gene))

blast_hits_rota <- blast_hits |> 
  group_by(sample_id) |>  
  summarise(unique_rota = paste(unique(gene), collapse = ", "))

write_xlsx(blast_hits_rota , "outer_capside_protein.xlsx")

######################################################################





