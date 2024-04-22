library(Biostrings)
library(XML)
library(seqinr)
library(dplyr)

#########################################################################
# Filter probes by blast
#########################################################################
# This script provides a comprehensive approach for managing probe data for BLAST analysis. 


# Function to load probes data, number them, and write to FASTA files
numbering_probes <- function(probes_rdata) {
  load(file = probes_rdata)
  if (exists("filtered_probes") && !is.null(filtered_probes)) {
    probes <- filtered_probes
  } else if (exists("fin_probes")) {
    probes <- fin_probes
  }
  probes_process <- probes
  probes_process <- probes %>%
    group_by(identifier) %>%
    mutate(probe_id = paste(gsub(" ", "_", identifier, fixed = TRUE), row_number(), sep = "_"))
  probes <- DNAStringSet(probes_process$probe[,1])
  target <- reverseComplement(DNAStringSet(probes))
  names(probes) <- probes_process$probe_id
  names(target) <- probes_process$probe_id
  writeXStringSet(probes, filepath = "./blast/probes.fasta")
  writeXStringSet(target, filepath = "./blast/target.fasta")
}

# Function to create BLAST database
create_blast_database <- function(makeblastdb_path,db_dir,db_name,query_db_path) {
  if (!dir.exists(db_dir)) {
    db_type <- "nucl"
    db_out <- paste(db_dir, db_name, sep = "/")
    dbmake <- paste(makeblastdb_path,
                    "-in", query_db_path,
                    "-dbtype", db_type,
                    "-out", db_out,
                    collapse = " ")
    system(dbmake)
    return(db_out)
  } else {
    cat("Database directory already exists. Skipping blast database creation.\n")
  }
}

# Function to run BLASTN
run_blastn <- function(blastn_path,query_blast,db,blast_out,word_size="5",evalue="1e-2",strand="plus",outfmt="5") {
  blast_cmd <- paste(blastn_path,
                     "-query",query_blast,
                     "-db", db,
                     "-word_size",word_size,
                     "-evalue",evalue,
                     "-strand",strand,
                     "-outfmt",outfmt,
                     "-out",blast_out,
                     collapse = " ")
  system(blast_cmd)
}

# Function to extract blstn xml file information 
xml_extract <- function(blast_xml_file){
  blast_xml <- xmlParse(blast_xml_file)
  hit_nodes <- getNodeSet(blast_xml, "//Iteration")
  blast_info_list <- lapply(hit_nodes, function(node) {
    query_id <- xmlValue(node[["Iteration_query-ID"]])
    hit_nodes <- getNodeSet(node, ".//Hit")
    hits <- lapply(hit_nodes, function(hit_node) {
      hit_id <- xmlValue(hit_node[["Hit_id"]])
      hit_title <- xmlValue(hit_node[["Hit_def"]])
      hsp_nodes <- getNodeSet(hit_node, ".//Hsp")
      hsps <- lapply(hsp_nodes, function(hsp_node) {
        hsp_align_len <- as.numeric(xmlValue(hsp_node[["Hsp_align-len"]]))
        hsp_identity <- as.numeric(xmlValue(hsp_node[["Hsp_identity"]]))
        hsp_gaps <- as.numeric(xmlValue(hsp_node[["Hsp_gaps"]]))
        hsp_query_seq <- xmlValue(hsp_node[["Hsp_qseq"]])
        hsp_hit_seq <- xmlValue(hsp_node[["Hsp_hseq"]])
        match_count <- hsp_identity
        matched_gc_count <- sum((unlist(strsplit(hsp_query_seq, "")) == unlist(strsplit(hsp_hit_seq, ""))) &
                                  unlist(strsplit(hsp_query_seq, "")) %in% c("G", "C"))
        return(data.frame(
          gap_count = hsp_gaps,
          match_count = match_count,
          matched_gc_count = matched_gc_count,
          stringsAsFactors = FALSE,
          hsp_query_seq = hsp_query_seq,
          hsp_hit_seq = hsp_hit_seq
        ))
      })
      match_count <- max(sapply(hsps,`[[`, "match_count"))
      matched_gc_count <- max(sapply(hsps,`[[`, "matched_gc_count"))
      hsp_query_seq = sapply(hsps,`[[`, "hsp_query_seq")
      hsp_hit_seq = sapply(hsps,`[[`, "hsp_hit_seq")
      return(data.frame(query_id = query_id,
                        hit_id = hit_id,
                        hit_title = hit_title,
                        match_count = match_count,
                        matched_gc_count = matched_gc_count,
                        stringsAsFactors = FALSE,
                        hsp_query_seq = hsp_query_seq,
                        hsp_hit_seq = hsp_hit_seq))
    })
    do.call(rbind, hits)
  })
  blast_info <- do.call(rbind, blast_info_list)
  return(blast_info)}

# Function to calculate eff between probes and hit_seqs
blast_calculate_eff <- function(blast_info, probes_fasta_path,temp=37,P=100e-9,ions=0.3,FA=50,batchSize=1000) {
  blast_info$probe_seq <- reverseComplement(DNAStringSet(blast_info$hsp_query_seq))

  probe_1 <- CalculateEfficiencyFISH(probe = blast_info$probe_seq,
                                     target = blast_info$hsp_hit_seq,
                                     temp,
                                     P,
                                     ions,   # Na
                                     FA,
                                     batchSize)
  selected_probe_1 <- probe_1[, 1:2]
  blast_eff <- cbind(blast_info, selected_probe_1)
  fasta_data <- readDNAStringSet(probes_fasta_path)
  fasta_df <- data.frame(
    probe_name = names(fasta_data),
    probe_seq_all = as.character(fasta_data),
    stringsAsFactors = FALSE
  )
  fasta_df$query_id <- paste("Query", 1:nrow(fasta_df), sep = "_")
  result_df <- merge(fasta_df, blast_eff, by = "query_id")
  result_df$probe_id <- gsub("_", " ", result_df$probe_name)

  probe_2 <- CalculateEfficiencyFISH(probe = result_df$probe_seq_all,
                                     target = result_df$hsp_hit_seq,
                                     temp,
                                     P,
                                     ions,   # Na
                                     FA,
                                     batchSize)
  selected_probe_2 <- probe_2[, 1:2]
  colnames(selected_probe_2) <- c("HybEff_alen", "FAm_alen")
  result_df <- cbind(result_df, selected_probe_2)
  
  return(result_df)
}

# Function to filter probes by match and eff
process_filtered_results <- function(result_df,max_match = 18,max_HybEff = 1,max_HybEff_alen = 0.01) {
  probe_match <- function(hit_title, probe_id) {
    target_species <- gsub(" \\d+$", "", probe_id)
    return(hit_title == target_species)
  }
  result_df$probe_seq <- as.character(result_df$probe_seq)

  filtered_df <- result_df %>%
    mutate(match = probe_match(hit_title, probe_id)) %>%
    group_by(probe_id) %>%
    filter(
      all((match == TRUE) | (match == FALSE   
                             & match_count <= max_match 
                             & HybEff < max_HybEff 
                             & HybEff_alen < max_HybEff_alen
      ))
    ) %>%
    ungroup()
  filtered_df$probe_id <- gsub(" \\d+$", "", filtered_df$probe_id)
  
  return(filtered_df)
}

# Function to merge all probes information
merged_results <- function(probes_rdata, filtered_df) {
  load(file = probes_rdata)
  if (exists("filtered_probes") && !is.null(filtered_probes)) {
    probes <- filtered_probes
  } else if (exists("fin_probes")) {
    probes <- fin_probes
  }
  probes <- probes %>%
    group_by(identifier) %>%
    mutate(probe_name = paste(gsub(" ", "_", identifier, fixed = TRUE), row_number(), sep = "_"))
  merge_df <- merge(filtered_df, probes, by = "probe_name")

  return(merge_df)
}

# Function to get probes subset by eff and aln position
group_results <- function(merged_df) {
  merge_df_sub <- merged_df %>%
    arrange(identifier, start) %>%
    group_by(identifier) %>%
    mutate(group = cumsum(c(0, diff(start) > 30))) %>%
    group_by(identifier, group) %>%
    slice(which.max(efficiency))
  return(merge_df_sub)
}

# Function to get probes subset with max_eff numbers
calculate_efficiency_stats <- function(grouped_df) {
  grouped_df <- grouped_df %>%
    group_by(identifier) %>%
    summarise(num_probes = n(),
              max_efficiency = max(efficiency, na.rm = TRUE),
              min_efficiency = min(efficiency, na.rm = TRUE),
              avg_efficiency = mean(efficiency, na.rm = TRUE)) %>%
    mutate(num_probes)
  
  return(grouped_df)
}

# Function to get two probes combinations if identifier only has one probe
process_sin_probe_id <- function(sub_probes_info, merge_df,gap_len=3) {
  sin_probe_ids <- sub_probes_info %>%
    filter(num_probes <= 1) %>%
    pull(identifier)
  merge_df <- merge_df %>% 
    arrange(identifier, start) %>%
    filter(identifier %in% sin_probe_ids)
  
  result_df <- data.frame(identifier = character(),
                          probe_name_1 = character(),
                          probe_name_2 = character(),
                          start_diff = numeric(),
                          start_1 = numeric(),
                          start_2 = numeric(),
                          avg_efficiency = numeric(),
                          HybEff_1 = numeric(),
                          HybEff_2 = numeric(),
                          probe_seq_1 = character(),
                          probe_seq_2 = character(),
                          stringsAsFactors = FALSE)
  
  for (id in unique(merge_df$identifier)) {
    id_probes <- merge_df %>% filter(identifier == id)
    num_probes <- nrow(id_probes)
    
    if (num_probes > 1) {
      probe_combinations <- list()
      for (i in 1:(num_probes - 1)) {
        for (j in (i + 1):num_probes) {
          probe1 <- id_probes[i, ]
          probe2 <- id_probes[j, ]

          if ((probe1$start + probe1$len + gap_len ) < probe2$start) {
            combination <- c(probe1$probe_name, probe2$probe_name)
            probe_combinations[[length(probe_combinations) + 1]] <- combination
          }
        }
      }
      
      if (length(probe_combinations) == 0) {
        cat("Identifier", id, "only has one probe.\n")
      } else {
        avg_efficiencies <- sapply(probe_combinations, function(comb) {
          probes <- merge_df %>% filter(probe_name %in% comb)
          mean(probes$HybEff)
        })
        
        for (k in 1:length(probe_combinations)) {
          if (length(probe_combinations[[k]]) >= 2) {
            combination <- probe_combinations[[k]]
            probe1 <- id_probes[id_probes$probe_name == combination[1], ]
            probe2 <- id_probes[id_probes$probe_name == combination[2], ]
            result_df <- rbind(result_df, data.frame(identifier = id,
                                                     probe_name_1 = combination[1],
                                                     probe_name_2 = combination[2],
                                                     start_diff = abs(probe1$start - probe2$start),
                                                     start_1 = probe1$start,
                                                     start_2 = probe2$start,
                                                     avg_efficiency = avg_efficiencies[k],
                                                     HybEff_1 = probe1$HybEff,
                                                     HybEff_2 = probe2$HybEff,
                                                     probe_seq_1 = probe1$probe,
                                                     probe_seq_2 = probe2$probe))
          }
        }
      }
    } else {
      cat("Identifier", id, "only has one probe.\n")
      result_df <- rbind(result_df, id_probes)
    }
  }
  
  return(result_df)
}

# Function to get best probes combination if identifier only has one probe
filter_best_combinations <- function(data,w1=0.5,w2=0.5){
  filtered_data <- data.frame(identifier = character(),
                              probe_name_1 = character(),
                              probe_name_2 = character(),
                              start_diff = numeric(),
                              avg_efficiency = numeric())
  identifiers <- unique(data$identifier)

  for (id in identifiers) {
    subset_data <- subset(data, identifier == id)
    weighted_score <- subset_data$start_diff * w1 + subset_data$avg_efficiency * w2 
    best_index <- which.max(weighted_score)
    best_combination <- subset_data[best_index, ]
    filtered_data <- rbind(filtered_data, best_combination)
  }
  
  return(filtered_data)
}

# Function to get and sort best two probes combination for all ids
select_best_probes_per_id <- function(merged_df, merge_df_sub, sin_best_combination) {
  multi_probe_ids <- merge_df_sub %>%
    group_by(identifier) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    select(identifier) %>%
    distinct()
  
  top_probes_multi <- merge_df_sub %>%
    filter(identifier %in% multi_probe_ids$identifier) %>%
    arrange(identifier, desc(HybEff)) %>%
    group_by(identifier) %>%
    slice_max(HybEff, n = 2, with_ties = FALSE) %>%
    ungroup()
  
  top_probes_multi <- top_probes_multi %>% select(-group)

  single_probe_ids <- merge_df_sub %>%
    group_by(identifier) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    select(identifier) %>%
    distinct() %>%
    pull(identifier)
  
  top_probes_single <- sin_best_combination %>%
    filter(identifier %in% single_probe_ids)
  
  get_probe_info <- function(probe_name, df) {
    df_row <- df[df$probe_name == probe_name, ]
    return(df_row)
  }
  
  top_probes_single_detailed <- data.frame()
  for (i in 1:nrow(top_probes_single)) {
    row <- top_probes_single[i, ]
    probe1_info <- get_probe_info(row$probe_name_1, merged_df)
    probe2_info <- get_probe_info(row$probe_name_2, merged_df)
    probe2_info <- probe2_info[colnames(probe1_info)]
    
    probe1_info$identifier <- row$identifier
    probe2_info$identifier <- row$identifier
    
    top_probes_single_detailed <- rbind(top_probes_single_detailed, probe1_info, probe2_info)
  }
  
  combined_probes_table <- bind_rows(top_probes_multi, top_probes_single_detailed)}

