set.seed(619291363)

library(magrittr)
library(DESeq2)

gene_counts_csv = "../data/swab_gene_counts.csv"
gene2name_txt = "../annotation/gene2name.txt"
cv_folds_csv = "../results/cv-folds.csv"

metatable_csv <- "../data/metatable_with_viral_status.csv"

## Read data

read.csv(metatable_csv, stringsAsFactors=F) %>%
  dplyr::inner_join(read.csv(cv_folds_csv, stringsAsFactors=F)) %>%
  dplyr::mutate(gender=make.names(gender)) %>%
  tibble::column_to_rownames("CZB_ID") %>%
  as.data.frame() ->
  metatable

gene_counts_csv %>%
  read.csv(row.names=1) %>%
  as.matrix() %>%
  .[,rownames(metatable)] ->
  gene_counts

gene2name_txt %>%
  read.csv(stringsAsFactors=F,
           col.names=c("ensembl_gene_id", "gene_name", "X"), sep="\t") %>%
  dplyr::select(ensembl_gene_id, gene_name) ->
  gene_names_df


## Outcome variable
y <- as.factor(dplyr::case_when(
  metatable$viral_status == "SC2" ~ "COVID19",
  metatable$viral_status == "other_virus" ~ "OtherVirus",
  metatable$viral_status == "no_virus" ~ "NoVirus",
  TRUE ~ NA_character_
))
covid <- metatable$viral_status == "SC2"

# folds
fold <- metatable$fold
uniq_folds <- unique(fold)

# training set indices for each fold
train_list <- lapply(uniq_folds, function(test_fold){ fold != test_fold })

# test y for each fold
test_y_list <- lapply(train_list, function(train) y[!train])

# VST

(metatable$gender == "M") %>%
  {. / sd(.)} ->
  gender_var

metatable$age %>%
  {. / sd(.)} ->
  age_var

X <- t(varianceStabilizingTransformation(gene_counts))
X_withAgeGender <- cbind(RenormalizedAge=age_var, RenormalizedGender=gender_var, X)

# Per-CV-fold VST

X_list_withAgeGender <- lapply(train_list, function(train) {
  dds_train <- DESeqDataSetFromMatrix(gene_counts[,train], metatable[train,], ~1)
  dds_train <- estimateSizeFactors(dds_train)
  dds_train <- estimateDispersions(dds_train)

  dds <- DESeqDataSetFromMatrix(gene_counts, metatable, ~1)
  dds <- estimateSizeFactors(dds)
  dispersionFunction(dds) <- dispersionFunction(dds_train)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

  X <- t(assay(vsd))
  stopifnot(rownames(X) == rownames(metatable))

  (metatable$gender == "M") %>%
    {. / sd(.[train])} ->
    gender_var

  metatable$age %>%
    {. / sd(.[train])} ->
    age_var

  cbind(RenormalizedAge=age_var, RenormalizedGender=gender_var, X)
})

X_list <- lapply(X_list_withAgeGender, function(X) X[,-(1:2)])


evaluate_predictions <- function(test_probs_list, test_y_list, out_prefix,
                                 exclude_category=NULL) {
  if (!is.null(exclude_category)) {
    test_probs_list <- mapply(function(probs, y) probs[y != exclude_category,],
                              test_probs_list, test_y_list, SIMPLIFY=FALSE)
    test_y_list <- lapply(test_y_list, function(y) y[y != exclude_category])
  }

  fold_rocs <- mapply(function(test_probs, test_y) {
    pROC::roc(test_y == "COVID19", test_probs[,"COVID19"])
  }, test_probs_list, test_y_list, SIMPLIFY=FALSE)

  sapply(fold_rocs, function(fr) fr$auc) %>%
    summary() %>%
    as.list() %>%
    cbind() %>%
    write.table(paste(out_prefix, "cv_auc_summary.txt", sep=""), col.names=F)

  sapply(fold_rocs, function(fr) fr$auc) %>%
    write(paste(out_prefix, "cv_auc.txt", sep=""), sep="\n")

  svg(paste(out_prefix, "cv_roc.svg", sep=""))
  for (i in 1:length(fold_rocs)) {
    plot(fold_rocs[[i]], add=(i != 1), col=i)
  }
  dev.off()

  test_preds_list <- lapply(test_probs_list, function(test_probs) {
    apply(test_probs, 1, function(row) names(row)[which.max(row)])
  })

  mapply(function(test_y, test_pred) {
    data.frame(truth=paste(test_y, ".true", sep=""),
               pred=paste(test_pred, ".pred", sep=""))
  }, test_y_list, test_preds_list, SIMPLIFY=FALSE) %>%
    do.call(what=rbind) %>%
    table() %>%
    write.csv(paste(out_prefix, "confusion_matrix.csv", sep=""))
}

evaluate_at_s <- function(X, X_list, mod, mod_list, s, out_dir) {
  dir.create(out_dir)

  probs_list <- mapply(function(glmnet_model, X, train){
    pred <- predict(glmnet_model, newx=X[!train,],
                    type="response", s=s)[,1]
    data.frame(COVID19= pred, NotCovid= 1 - pred)
  }, mod_list, X_list, train_list, SIMPLIFY=FALSE)

  evaluate_predictions(
    probs_list, test_y_list, file.path(out_dir, "lasso_"))

  evaluate_predictions(
    probs_list, test_y_list,
    file.path(out_dir, "lasso_excludeOther_"),
    exclude_category="OtherVirus")

  evaluate_predictions(
    probs_list, test_y_list,
    file.path(out_dir, "lasso_viralOnly_"),
    exclude_category="NoVirus")

  coef(mod, s=s)[,1] %>%
    .[. != 0] ->
    nonzero_coefs

  nonzero_coefs %>%
    {data.frame(variable=names(.), coef=.)} %>%
    {dplyr::right_join(
      dplyr::rename(gene_names_df, variable=ensembl_gene_id), .)} %>%
    write.csv(file.path(out_dir, "lasso_nonzero_coefs.csv"), row.names=F)

  X_subsetColumns_list <- mapply(function(glmnet_model, X) {
    X[, names(nonzero_coefs)[-1]]
  }, mod_list, X_list, SIMPLIFY=FALSE)

  rf_list <- mapply(function(X, train) {
    randomForest::randomForest(
      X[train,], as.factor(dplyr::if_else(covid[train], "COVID19", "NotCovid")), ntree=10000)
  }, X_subsetColumns_list, train_list, SIMPLIFY=FALSE)

  rf_probs_list <- mapply(function(rf, X, train){
    predict(rf, newdata=X[!train,], type="prob")
  }, rf_list, X_subsetColumns_list, train_list, SIMPLIFY=FALSE)

  mapply(function(rf_probs, train){
    data.frame(CZB_ID=row.names(metatable)[!train], rf_probs,
               stringsAsFactors=FALSE)
  }, rf_probs_list, train_list, SIMPLIFY=FALSE) %>%
    do.call(what=rbind) %>%
    write.csv(file.path(out_dir, "lassoRandomForest_probs.csv"), row.names=FALSE)

  evaluate_predictions(
    rf_probs_list, test_y_list,
    file.path(out_dir, "lassoRandomForest_"))

  evaluate_predictions(
    rf_probs_list, test_y_list,
    file.path(out_dir, "lassoRandomForest_excludeOther_"),
    exclude_category="OtherVirus")

  evaluate_predictions(
    rf_probs_list, test_y_list,
    file.path(out_dir, "lassoRandomForest_viralOnly_"),
    exclude_category="NoVirus")
}

out_dir <- "../results/classifier"
dir.create(out_dir)

lasso_model <- function(X, train, penalty.factor) {
  glmnet::cv.glmnet(
    X[train,], covid[train],
    foldid=as.integer(as.factor(fold[train])),
    penalty.factor=penalty.factor,
    standardize=FALSE, family="binomial")
}

covid_mod_list <- mapply(
  lasso_model, X_list, train_list,
  MoreArgs=list(penalty.factor=rep(1, ncol(X))), SIMPLIFY=FALSE)

covid_mod <- lasso_model(X, rep(TRUE, nrow(X)), rep(1, ncol(X)))

withAgeGender_mod_list <- mapply(
  lasso_model, X_list_withAgeGender, train_list,
  MoreArgs=list(penalty.factor=c(0, 0, rep(1, ncol(X)))), SIMPLIFY=FALSE)

withAgeGender_mod <- lasso_model(X_withAgeGender, rep(TRUE, nrow(X)), c(0, 0, rep(1, ncol(X))))

evaluate_at_s(X, X_list, covid_mod, covid_mod_list, "lambda.min", file.path(out_dir, "lasso_min"))
evaluate_at_s(X, X_list, covid_mod, covid_mod_list, "lambda.1se", file.path(out_dir, "lasso_1se"))

s10 <- min(covid_mod$lambda[covid_mod$nzero <= 10])
evaluate_at_s(X, X_list, covid_mod, covid_mod_list, s10, file.path(out_dir, "lasso_s10"))

s5 <- min(covid_mod$lambda[covid_mod$nzero <= 5])
evaluate_at_s(X, X_list, covid_mod, covid_mod_list, s5, file.path(out_dir, "lasso_s5"))

s3 <- min(covid_mod$lambda[covid_mod$nzero <= 3])
evaluate_at_s(X, X_list, covid_mod, covid_mod_list, s3, file.path(out_dir, "lasso_s3"))

evaluate_at_s(X_withAgeGender, X_list_withAgeGender, withAgeGender_mod, withAgeGender_mod_list, "lambda.min", file.path(out_dir, "withAgeGender_lasso_min"))
evaluate_at_s(X_withAgeGender, X_list_withAgeGender, withAgeGender_mod, withAgeGender_mod_list, "lambda.1se", file.path(out_dir, "withAgeGender_lasso_1se"))

s10 <- min(withAgeGender_mod$lambda[withAgeGender_mod$nzero <= 10])
evaluate_at_s(X_withAgeGender, X_list_withAgeGender, withAgeGender_mod, withAgeGender_mod_list, s10, file.path(out_dir, "withAgeGender_lasso_s10"))

s5 <- min(withAgeGender_mod$lambda[withAgeGender_mod$nzero <= 5])
evaluate_at_s(X_withAgeGender, X_list_withAgeGender, withAgeGender_mod, withAgeGender_mod_list, s5, file.path(out_dir, "withAgeGender_lasso_s5"))

s3 <- min(withAgeGender_mod$lambda[withAgeGender_mod$nzero <= 3])
evaluate_at_s(X_withAgeGender, X_list_withAgeGender, withAgeGender_mod, withAgeGender_mod_list, s3, file.path(out_dir, "withAgeGender_lasso_s3"))
