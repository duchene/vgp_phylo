#!/usr/bin/env Rscript
## Stage 0 covariate: the chromosome a locus lies on in two reference species, from
## the sequence labels (results/stage0/locus_labels.tsv.gz) and the NCBI assembly
## reports (results/stage0/assembly_reports/). Membership only: the labels give the
## chromosome, not the position on it.
##
##   chicken (GCF_016700215.2): chromosome, its length, and a size class
##     macro = 1-5, Z, W; intermediate = 6-10; micro = 11 and up (docs/PLAN, Step 1)
##   human (GCA_000001405.15): chromosome and its length (the deep clades)
##
## Loci without the species, or whose sequence sits on an unplaced scaffold, get NA.
## Usage: Rscript scripts/stage0_chromosomes.R -> results/stage0/locus_chromosome.tsv

root <- Sys.getenv("VGP_PHYLO_ROOT", unset = normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
source(file.path(root, "config", "project.R"))
S0 <- file.path(PATHS$root, "results", "stage0")

REF <- c(chicken = "GCF_016700215.2", human = "GCA_000001405.15")
lab <- utils::read.delim(gzfile(file.path(S0, "locus_labels.tsv.gz")), stringsAsFactors = FALSE)
lab <- lab[lab$assembly %in% REF, ]
m <- utils::read.delim(PATHS$manifest, stringsAsFactors = FALSE)
out <- data.frame(locus_id = m$locus_id, clade = m$clade, stringsAsFactors = FALSE)

for (sp in names(REF)) {
  rep <- utils::read.delim(file.path(S0, "assembly_reports", paste0(REF[[sp]], ".tsv")),
                           stringsAsFactors = FALSE, colClasses = "character")
  rep <- rep[rep$assigned_molecule_location_type == "Chromosome" & rep$role == "assembled-molecule", ]
  ## labels carry the GenBank or the RefSeq accession, or (human) a UCSC name "chr2"
  key <- c(rep$genbank_accession, rep$refseq_accession, paste0("chr", rep$chr_name))
  chr <- rep(rep$chr_name, 3); len <- rep(as.numeric(rep$length), 3)
  l <- lab[lab$assembly == REF[[sp]], ]
  if (anyDuplicated(l$locus_id)) stop(sp, ": more than one sequence per locus")
  i <- match(out$locus_id, l$locus_id); j <- match(l$seq_accession[i], key)
  out[[paste0(sp, "_chr")]] <- chr[j]
  out[[paste0(sp, "_chr_len")]] <- len[j]
  cat(sprintf("%-8s present in %d loci, on an assembled chromosome in %d\n", sp,
              sum(!is.na(i)), sum(!is.na(j))))
}
out$chicken_class <- with(out, ifelse(is.na(chicken_chr), NA,
  ifelse(chicken_chr %in% c(as.character(1:5), "Z", "W"), "macro",
  ifelse(chicken_chr %in% as.character(6:10), "intermediate", "micro"))))
print(table(out$chicken_class[out$clade == "BirdsAnc0"], useNA = "ifany"))

f <- file.path(S0, "locus_chromosome.tsv")
utils::write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(nrow(out), "loci ->", f, "\n")
