library(data.table)
library(magrittr)

Template_alleles <- fread("Template_DB_information.csv")[, MHC] %>% unique

freq_table <- fread("./mhcflurry.ba.frequency_matrices.csv")
freq_table <- freq_table[cutoff_fraction == 0.001, !c("cutoff_fraction", "cutoff_count")]

HELL <- function(c1, c2) {
  1 - mean(sqrt(rowSums((sqrt(c1) - sqrt(c2))^2)) / sqrt(2))
}

extract_matrix <- function(freq_table, allele_list, Template_alleles, nmer) {
  allele_list <- split(freq_table, by="allele", keep.by=FALSE)
  allele_list <-lapply(allele_list, function(mat) mat[,-1])
  alleles_found_in_MHCFlurry <- names(allele_list)[names(allele_list) %in% Template_alleles]
  not_supported_alleles <- setdiff(Template_alleles, names(allele_list)[names(allele_list) %in% Template_alleles])
  filtered_allele_list <- allele_list[alleles_found_in_MHCFlurry]
  i <- 0
  similarity_matrix <- lapply(allele_list, function(x) {
    similarity <- sapply(filtered_allele_list, function(y) {HELL(x, y)})
    similarity <- c(similarity, rep(0, length(not_supported_alleles))) # Min possible Hellinger similarity 
    i <<- i + 1
    message(paste(i, "/", length(allele_list)),"\r",appendLF=FALSE)
    flush.console()
    similarity
  })
  similarity_matrix <- do.call(rbind, similarity_matrix)
  colnames(similarity_matrix) <- c(names(filtered_allele_list), not_supported_alleles)
  similarity_matrix <- rbind(similarity_matrix, matrix(0, nrow = length(not_supported_alleles), ncol = length(Template_alleles)))
  rownames(similarity_matrix) <- c(names(allele_list), not_supported_alleles)
  similarity_matrix[substr(rownames(similarity_matrix), 1, 4) == "Gaga", substr(colnames(similarity_matrix), 1, 4) == "Gaga"] <- 0.5
  similarity_matrix[outer(rownames(similarity_matrix), colnames(similarity_matrix), "==")] <- 1 # Max possible Hellinger similarity 
  similarity_df <- as.data.table(similarity_matrix, keep.rownames = TRUE)
  setnames(similarity_df, "rn", "Allele")
  fwrite(similarity_df, paste0(nmer, "mer_similarity.csv"))
}

for (nmer in seq(8,15)) {
  extract_matrix(freq_table[length == nmer, !c("length")], allele_list, Template_alleles, nmer)
}

#### CAUTION!!! Go and increase similarity between chciken alleles manually after this!

# What I want to do here is calculate 

# Let's see if I can make a tree:
#correlation_matrix[alleles_found_in_MHCFlurry, alleles_found_in_MHCFlurry]

#hc <- hclust(as.dist(correlation_matrix[alleles_found_in_MHCFlurry, alleles_found_in_MHCFlurry]), method = "ward.D2")
#plot(unroot(as.phylo(hc)),type="unrooted",cex=1,use.edge.length=FALSE,lab4ut="axial",no.margin=TRUE,direction="downwards", rotate.tree = 0, main="Scaled Phylogenetic tree")
