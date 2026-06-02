# Merge count data generated from BAM file with the blast data
# Load packages:
library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(writexl)
library(tidyverse)

# Files to combine
virus_file  <- "ncbi_dataset/data/rota_virus.tsv" # This file has al the description
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

write_xlsx(merged, "rota_count_blast.xlsx")

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
