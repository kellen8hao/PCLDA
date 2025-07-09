library(MASS)  # for LDA
library(caret) # for confusionMatrix (optional)


set.seed(1)

commandIn <- commandArgs(F)
print(commandIn)
cmdArgs <- commandArgs(TRUE) ## XJ: reads command from the shell


train <- readRDS(cmdArgs[1]); 
test <- readRDS(cmdArgs[2]);
label_output_file <- cmdArgs[3]; #"pre_lable" 
aa_output_file <- cmdArgs[4]; # "accuracy"


cell_label=train[,1]
test_y=test[,1]


### step 2: build result matrix
a_mtx <- matrix(NA, 1, 3)
colnames(a_mtx)=c("dataset","accuracy","time")
a_mtx[,1]=basename(dirname(dirname(cmdArgs[1])))


label_mtx<- matrix(NA, length(test_y), 2)
rownames(label_mtx)=rownames(test)
colnames(label_mtx)=c("yy_int","yy_pred")





run_pclda <- function(train_matrix, test_matrix, cell_label, test_y, num_of_pc = 200, scale_data = TRUE) {
  start=Sys.time()
  # Step 1: Optionally scale data
  if (scale_data) {
    train_scaled <- scale(train_data[,-1])
    test_scaled <- scale(test_data[,-1], center = attr(train_scaled, "scaled:center"), scale = attr(train_scaled, "scaled:scale"))
  } else {
    train_scaled <- train_data[,-1]
    test_scaled <- test_data[,-1]
  }
  
  # Step 1: PCA on combined matrix
  new_matrix = rbind(train_matrix, test_matrix)
  pca_re = prcomp(new_matrix, scale = FALSE)
  eigvecs = -1 * pca_re$rotation
  
  # Step 2: Compute total variance (between-class)
  total_class_vars = var(train_matrix)
  
  # Step 3: Compute within-class variance
  celltypes = unique(cell_label)
  within_class_vars = var(train_matrix[cell_label == celltypes[[1]], , drop = FALSE])
  for (ct in celltypes[-1]) {
    within_class_vars = within_class_vars + var(train_matrix[cell_label == ct, , drop = FALSE])
  }
  
  # Step 4: Compute P_within and P_total for each PC
  p_within_metric <- function(eigvec) {
    t(eigvec) %*% within_class_vars %*% eigvec
  }
  p_total_metric <- function(eigvec) {
    t(eigvec) %*% total_class_vars %*% eigvec
  }
  p_within = apply(eigvecs, 2, p_within_metric)
  p_total  = apply(eigvecs, 2, p_total_metric)
  p_ratio  = p_total / p_within
  
  # Step 5: Select top PCs by sorted ratio
  sorted = sort(p_ratio, index.return = TRUE, decreasing = TRUE)
  selected_eigvecs = eigvecs[, sorted$ix[1:num_of_pc]]
  
  # Step 6: Project train/test onto selected PCs
  projected_expr_matrix = as.matrix(train_matrix) %*% selected_eigvecs * -1
  projected_test_matrix = as.matrix(test_matrix) %*% selected_eigvecs * -1
  
  # Step 7: Train LDA
  train_pca = data.frame(yy = cell_label, projected_expr_matrix)
  test_pca  = data.frame(yy = test_y, projected_test_matrix)
  formula_str = paste("yy ~", paste(colnames(projected_expr_matrix), collapse = "+"))
  lda_fit = lda(as.formula(formula_str), data = train_pca, tol = 1e-20)
  
  # Step 8: Predict
  lda_pred = predict(lda_fit, test_pca)
  predicted_labels = lda_pred$class
  
  end = Sys.time()
  
  # Output
  return(list(
    predicted = predicted_labels,
    posterior = lda_pred$posterior,
    model = lda_fit,
    time = difftime(end, start, units = "secs"),
    p_ratio = p_ratio,
    selected_pcs = colnames(projected_expr_matrix)
  ))
}


## output

result <- pclda_classifier(train_matrix, test_matrix, cell_label, test_y, num_of_pc = 200, scale_data = TRUE)
head(result$predicted)

## accuracy
pre_label=result$predicted
lab_c=table(pre_label==test_y)["TRUE"] 
accuracy=lab_c/length(test_y)
accuracy
print("accuary for lda  with ratio:")

## 5. save result
label_mtx[,1]=test_y
label_mtx[,2] <- as.character(pre_label)

a_mtx[,2]=accuracy
a_mtx[,3]=result$time


saveRDS(label_mtx, label_output_file)
write.csv(a_mtx, aa_output_file)

