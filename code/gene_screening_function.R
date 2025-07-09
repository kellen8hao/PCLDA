library(dplyr)
library(Matrix)
library(MASS)

commandIn <- commandArgs(F)
print(commandIn)
cmdArgs <- commandArgs(TRUE) ## XJ: reads command from the shell

train <- readRDS(cmdArgs[1]); 
test <- readRDS(cmdArgs[2]); 
train_output_file <- cmdArgs[3]; 
test_output_file <- cmdArgs[4]; 



select_top_genes_by_ttest <- function(train_data, test_data, top_n = 400, output_train_path = NULL, output_test_path = NULL) {
  start_time <- Sys.time()
  
  # Step 1: Match common genes
  common_genes <- intersect(colnames(train_data), colnames(test_data))
  train_data <- train_data[, common_genes]
  test_data <- test_data[, common_genes]
  
  celltypes <- unique(train_data[, 1])
  gene_names <- colnames(train_data)[-1]
  num_genes <- length(gene_names)
  num_ct <- length(celltypes)
  
  # Initialize p-value matrix
  pval_matrix <- matrix(NA, nrow = num_genes, ncol = num_ct)
  rownames(pval_matrix) <- gene_names
  colnames(pval_matrix) <- celltypes
  
  for (ct in seq_along(celltypes)) {
    is_ct <- train_data[, 1] == celltypes[ct]
    for (g in seq_along(gene_names)) {
      x <- train_data[is_ct, g + 1]
      y <- train_data[!is_ct, g + 1]
      t_res <- tryCatch(t.test(x, y, paired = FALSE), error = function(e) list(p.value = 1))
      pval_matrix[g, ct] <- ifelse(is.na(t_res$p.value), 1, t_res$p.value)
    }
  }
  
  # Step 2: Select top N genes per cell type
  selected_genes <- unique(unlist(apply(pval_matrix, 2, function(p) names(sort(p)[1:top_n]))))
  
  # Step 3: Subset and match again
  selected_cols <- c("yy", selected_genes)
  train_filtered <- data.frame(yy = train_data[, 1], train_data[, selected_genes, drop = FALSE])
  test_filtered <- data.frame(yy = test_data[, 1], test_data[, selected_genes, drop = FALSE])
  
  # Ensure same column order
  common_final <- intersect(colnames(train_filtered), colnames(test_filtered))
  train_filtered <- train_filtered[, common_final]
  test_filtered <- test_filtered[, common_final]
  
  end_time <- Sys.time()
  cat("Time for t-test selection:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  cat("Train dim after filtering:", dim(train_filtered), "\n")
  cat("Test dim after filtering:", dim(test_filtered), "\n")
  
  # Optional save
  if (!is.null(output_train_path)) saveRDS(train_filtered, output_train_path)
  if (!is.null(output_test_path)) saveRDS(test_filtered, output_test_path)
  
  return(list(train = train_filtered, test = test_filtered, selected_genes = selected_genes))
}


result <- select_top_genes_by_ttest(
  train_data = train,
  test_data = test,
  top_n = 400,
  output_train_path = train_output_file,
  output_test_path = test_output_file
)

# Access results
filtered_train <- result$train
filtered_test <- result$test
selected_genes <- result$selected_genes

