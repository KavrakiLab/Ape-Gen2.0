library(data.table)
library(magrittr)
#library(ade4)
#library(ape)
#library(dendextend)
#library(ggplot2)
#library(ggrepel)

Template_alleles <- fread("Template_information.csv")[, MHC] %>% unique

freq_table <- fread("./mhcflurry.ba.frequency_matrices.csv")
freq_table <- freq_table[cutoff_fraction == 0.001, !c("cutoff_fraction", "cutoff_count")]

freq_table_8 <- freq_table[length == 8, !c("length")]
freq_table_9 <- freq_table[length == 9, !c("length")]
freq_table_10 <- freq_table[length == 10, !c("length")]
freq_table_11 <- freq_table[length == 11, !c("length")]
freq_table_12 <- freq_table[length == 12, !c("length")]
freq_table_13 <- freq_table[length == 13, !c("length")]
freq_table_14 <- freq_table[length == 14, !c("length")]
freq_table_15 <- freq_table[length == 15, !c("length")]

extract_matrix <- function(freq_table, allele_list, Template_alleles, nmer) {
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  allele_list <- split(freq_table, by="allele", keep.by=FALSE)
  alleles_found_in_MHCFlurry <- names(allele_list)[names(allele_list) %in% Template_alleles]
  not_supported_alleles <- setdiff(Template_alleles, names(allele_list)[names(allele_list) %in% Template_alleles])
  filtered_allele_list <- allele_list[alleles_found_in_MHCFlurry]
  i <- 0
  correlation_matrix <- lapply(allele_list, function(x) {
    correlation <- sapply(filtered_allele_list, function(y) {
      euclidean(x %>% as.matrix %>% t %>% as.vector, y %>% as.matrix %>% t %>% as.vector)
    })
    correlation <- c(correlation, rep(180, length(not_supported_alleles))) # Max possible euclidean distance 
    i <<- i + 1
    message(paste(i, "/", length(allele_list)),"\r",appendLF=FALSE)
    flush.console()
    correlation
  })
  correlation_matrix <- do.call(rbind, correlation_matrix)
  colnames(correlation_matrix) <- c(names(filtered_allele_list), not_supported_alleles)
  correlation_matrix <- rbind(correlation_matrix, matrix(180, nrow = length(not_supported_alleles), ncol = length(Template_alleles)))
  rownames(correlation_matrix) <- c(names(allele_list), not_supported_alleles)
  correlation_matrix[outer(rownames(correlation_matrix), colnames(correlation_matrix), "==")] <- 0
  correlation_df <- as.data.table(correlation_matrix, keep.rownames = TRUE)
  setnames(correlation_df, "rn", "Allele")
  fwrite(correlation_df, paste0(nmer, "mer_similarity.csv"))
}

extract_matrix(freq_table_8, allele_list, Template_alleles, "8")
extract_matrix(freq_table_9, allele_list, Template_alleles, "9")
extract_matrix(freq_table_10, allele_list, Template_alleles, "10")
extract_matrix(freq_table_11, allele_list, Template_alleles, "11")
extract_matrix(freq_table_12, allele_list, Template_alleles, "12")
extract_matrix(freq_table_13, allele_list, Template_alleles, "13")
extract_matrix(freq_table_14, allele_list, Template_alleles, "14")
extract_matrix(freq_table_15, allele_list, Template_alleles, "15")


# Let's see if I can make a tree:
#correlation_matrix[alleles_found_in_MHCFlurry, alleles_found_in_MHCFlurry]

#hc <- hclust(as.dist(correlation_matrix[alleles_found_in_MHCFlurry, alleles_found_in_MHCFlurry]), method = "ward.D2")
#plot(unroot(as.phylo(hc)),type="unrooted",cex=1,use.edge.length=FALSE,lab4ut="axial",no.margin=TRUE,direction="downwards", rotate.tree = 0, main="Scaled Phylogenetic tree")
