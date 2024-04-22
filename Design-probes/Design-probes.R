library(DBI)
library(tools)
library(RSQLite)
library(Biostrings)
library(DECIPHER)

#########################################################################
# Design probes
#########################################################################
# This script processes 16S rRNA sequences for probe design and validation.
# It integrates sequence manipulation, taxonomic identification, and probe specificity testing against a reference database.
# Utilizes 'DBI' for database operations.
# Uses 'Biostrings' and 'DECIPHER' for sequence analysis and probe design.
# Employs 'tools' for file management.
# Ensures efficient and accurate probe development for microbial identification.

# Function to process raw datas
raw_16S_seq_process <- function(folder_path){
  # fasta_file <- file.path(folder_path, "16S.fasta")
  # if (file.exists(fasta_file)) {
  #   message("16S.fasta already exists in the specified folder. Skipping processing.")
  #   return()
  # }
  files <- list.files(folder_path, full.names = TRUE) 
  # suffix_to_remove <- "-final"
  rRNA_seqs <- vector("list", length = length(files))
  names_list <- vector("character", length = length(files))
  for (i in seq_along(files)) {
    file <- files[i]
    sequence <- paste0(readLines(file, warn = FALSE), collapse = "")
    base_name <- file_path_sans_ext(basename(file))
    #name <- gsub(suffix_to_remove, "", base_name)
    rRNA_seqs[[i]] <- sequence
    #names_list[i] <- name
    names_list[i] <- base_name
  }
  rRNA_seqs <- DNAStringSet(unlist(rRNA_seqs))
  names(rRNA_seqs) <- names_list
  writeXStringSet(rRNA_seqs,filepath = "./16S.fasta")}

# Function to get all bac taxonomy info
get_taxonomy_info <- function(input_fasta,taxonkit_path="/home/my8rose/anaconda3/bin/taxonkit") {
  non_standard_names <- FALSE  
  if (!file.exists("name_taxid.txt")) {
    if (!file.exists("name.txt")) {
      seqs <- readDNAStringSet(filepath = input_fasta)
      seq_names <- names(seqs)
      writeLines(seq_names, "name.txt")
    }
    
    input_data <- readLines("name.txt")
    output_data <- system2(taxonkit_path, args = "name2taxid", input = input_data, stdout = TRUE)
    for (i in seq_along(output_data)) {
      split_line <- strsplit(output_data[i], "\t")[[1]]
      if (length(split_line) < 2) {
        non_standard_names <- TRUE
        cat(input_data[i], "naming format is not standardized\n")
      }
    }
    if (non_standard_names) {
      return()
    }
    writeLines(output_data, "name_taxid.txt")
  }
  
  name_taxid <- read.table("name_taxid.txt", sep = "\t", stringsAsFactors = FALSE, quote = "")
  taxid <- name_taxid[,2]
  write.table(taxid, "taxid.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  input_data <- readLines("taxid.txt")
  output_data <- system2(taxonkit_path, args = "lineage", input = input_data, stdout = TRUE)
  writeLines(output_data, "lineage.txt")
  
  lineage <- readLines("lineage.txt")
  system2(taxonkit_path, args = c("reformat"), 
          input = lineage, stdout = "lineage_temp.txt")
  lineagetemp <- read.table(file = "lineage_temp.txt",sep="\t")
  lineagetemp <- lineagetemp[,-2]
  colnames(lineagetemp)[1] <- "taxid"
  split_columns <- strsplit(lineagetemp$V3, ";")
  max_columns <- max(sapply(split_columns, length))
  lineagetemp_split <- cbind(lineagetemp, t(sapply(split_columns, function(x) {
    length_diff <- max_columns - length(x)
    if(length_diff > 0) {
      x <- c(x, rep(NA, length_diff))
    }
    return(x)
  })))
  lineage_split<- lineagetemp_split[,-2]
  colnames(lineage_split)[2:8]<- c("kingdom","phylum","class","order","family","genus","species")
  colnames(name_taxid) <- c("identifier","rtaxid")
  merged_data <- merge(name_taxid, lineage_split, by.x = "rtaxid", by.y = "taxid", all.x = TRUE)
  
  write.table(merged_data, "16S_taxon_info.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  taxon_info <- read.table(file = "./16S_taxon_info.txt",sep = "\t",header = TRUE)
  
  file.remove("taxid.txt", "lineage.txt", "lineage_temp.txt")
  return(taxon_info)
}

# Function to design probes
process_probes_design_own <- function(input_fasta,
                                      db_file,
                                      minLength_tiles = 28,
                                      maxLength_tiles = 40,
                                      maxTilePermutations = 10,
                                      minLength_dp=28,
                                      maxLength_dp=34,
                                      maxPermutations = 4,
                                      hybTemp = 37,
                                      P = 100e-9,
                                      Na = 0.3,
                                      FA = 50){
  
  rRNA_seqs <- readDNAStringSet(filepath = input_fasta)
  con <- dbConnect(SQLite(), dbname = db_file)
  Seqs2DB(seqs = rRNA_seqs,
          type = "XStringSet",
          dbFile = con,
          identifier = names(rRNA_seqs),
          tblName = "Seqs",
          verbose = TRUE,
          replaceTbl = FALSE,
          processors = 8)
  tiles <- TileSeqs(con, 
                    add2tbl="Tiles",
                    identifier = "",
                    minLength = minLength_tiles,
                    maxLength = maxLength_tiles,
                    processors = 10)
  dbDisconnect(con)
  probes <- DesignProbes(tiles,
                         identifier = "",
                         minLength = minLength_dp,
                         maxLength = maxLength_dp,
                         maxPermutations = maxPermutations,
                         hybTemp = hybTemp,
                         P = P,
                         Na = Na, 
                         FA = FA,
                         verbose = TRUE)
  #probes <- probes %>% select_if(~ !any(is.na(.)))
  return(probes)}

# Function to filter probes by GC and dG
process_probes_filter <- function(probes,
                                  minLength_dp = 28,
                                  maxLength_dp = 34,
                                  min_GC = 0.45,
                                  max_GC = 0.65,
                                  temp = 37,
                                  P = 100e-9,
                                  ions = 0.3,
                                  FA = 50,
                                  batchSize = 1000,
                                  max_dG = -13.0){
  
  probes$len <- nchar(probes$probe[,1]) 
  probes <- probes[probes$len >= minLength_dp & probes$len <= maxLength_dp,]
  probes$GC_content <- sapply(probes$probe[,1], function(seq) {
    gc_count <- sum(sapply(c("G", "C"), function(base) {
      length(gregexpr(pattern = base, text = seq, ignore.case = TRUE)[[1]])
    }))
    gc_content <- gc_count / nchar(seq)
    return(gc_content)
  })
  
  filtered_probes <- probes[probes$GC_content >= min_GC & probes$GC_content <= max_GC, ]
  probe <- filtered_probes$probe[,1]
  target <- reverseComplement(DNAStringSet(probe))
  probe_G <- CalculateEfficiencyFISH(probe,
                                     target,
                                     temp,
                                     P,
                                     ions,   #Na
                                     FA,
                                     batchSize)
  # head(probe_G)
  filtered_probes$dG1 <- probe_G[,4]
  filtered_probes <- filtered_probes[filtered_probes$dG1 < max_dG, ]
  return(filtered_probes)}

# Function to get ref database tiles
ref_database_info <- function(ref_db_file = "~/R_data/ref_16S.sqlite"){
  dbConn_ref <- dbConnect(SQLite(), ref_db_file)
  ref_tiles <- dbGetQuery(dbConn_ref,
                          "select * from Tiles where groupCoverage > 0.2 and coverage > 0.01")
  dbDisconnect(dbConn_ref)
  ref_tiles$id <- paste("ref", ref_tiles$id, sep="_")
  return(ref_tiles)}

# Function to align the probes with the ref database tiles
compare_probes_refdatabase <- function(ref_tiles,
                                       probes,
                                       Hyb_FA=50,
                                       hybTemp=37,
                                       P_con=100e-9,
                                       Na_con=0.3,
                                       taxon_info){
  # The probes are compared to a reference database
  seqs <- DNAStringSet(ref_tiles$target_site)
  w <- which(!is.na(t(probes$probe)))
  probes_rc <- reverseComplement(DNAStringSet(t(probes$probe)[w]))
  p <- PDict(probes_rc, tb.start=1, tb.width=5)
  hits1 <- vwhichPDict(p, seqs, max.mismatch=5)
  l <- vapply(hits1, length, integer(1))
  hits1 <- unlist(hits1, use.names=FALSE)
  names(hits1) <- rep(1:length(l), l)
  p <- PDict(probes_rc, tb.end=-1, tb.width=5)
  hits2 <- vwhichPDict(p, seqs, max.mismatch=5)
  l <- vapply(hits2, length, integer(1))
  hits2 <- unlist(hits2, use.names=FALSE)
  names(hits2) <- rep(1:length(l), l)
  hits <- c(hits1, hits2)
  
  #Hyb_FA <- 50 
  count <- 0L 
  pBar <- txtProgressBar(style = 3) 
  for (i in 1:dim(probes)[1]) {
    # Extract the genus corresponding to the current probe
    current_genus <- taxon_info$genus[which(taxon_info$identifier == probes$identifier[i])]
    w <- which(!is.na(probes[i, "probe"]))
    results <- NULL
    for (j in w) {
      count <- count + 1L
      w_hits <- which(hits == count)
      if (length(w_hits) > 0) {
        ns <- as.integer(unique(names(hits[w_hits])))
        valid_ns <- ns 
        
        # Filter out reference sequences that are the same as the current probe genus
        for (n in ns) {
          ref_genus <- gsub("^ref_([^_]+)_?\\d*$", "\\1", ref_tiles$id[n])
          if (ref_genus == current_genus) {
            valid_ns <- valid_ns[valid_ns != n]
          }
        }
        
        if (length(valid_ns) > 0) {
          eff <- CalculateEfficiencyFISH(rep(probes$probe[i, j], length(valid_ns)),
                                         ref_tiles$target_site[valid_ns],
                                         hybTemp, 
                                         P_con, 
                                         Na_con, 
                                         Hyb_FA)
          
          eff <- cbind(eff,
                       data.frame(id=ref_tiles$id[valid_ns],
                                  probe_rc=toString(probes_rc[count]),
                                  target=ref_tiles$target_site[valid_ns],
                                  dFAm=eff[, "FAm"] - Hyb_FA,
                                  stringsAsFactors=FALSE))
          results <- rbind(results, eff)
        }
      }
    }
    w <- which(results$dFAm > -20)
    if (length(w) > 0) {
      # only record the strongest cross-hybridization in each group
      results <- results[w,]
      u <- unique(results$id)
      keep <- integer()
      for (j in 1:length(u)) {
        w <- which(results$id==u[j])
        keep <- c(keep, w[which.max(results$dFAm[w])])
      }
      results <- results[keep,]
      
      # append more non-targets to mismatches
      p <- pairwiseAlignment(results$probe_rc,
                             results$target,
                             type="global-local")
      probes$mismatches[i] <- paste(probes$mismatches[i],
                                    paste(results$id,
                                          " (", round(100*results$HybEff, 1), "%,",
                                          round(results$ddG1, 2), "kcal/mol,",
                                          round(results$dFAm, 1), "%;",
                                          substring(reverseComplement(DNAStringSet(pattern(p))),
                                                    1L),
                                          "/", substring(subject(p), 1L), ")",
                                          sep="",
                                          collapse=" "),
                                    sep="")
      # score -= 0.2 + 1.2^dFAm
      probes$score[i] <- probes$score[i] -
        sum(ifelse(results$dFAm < -20, 0, 0.2 +
                     1.2^ifelse(results$dFAm > 0, 0, results$dFAm)))
    }
    setTxtProgressBar(pBar, i/nrow(probes))
  }
  return(probes)}