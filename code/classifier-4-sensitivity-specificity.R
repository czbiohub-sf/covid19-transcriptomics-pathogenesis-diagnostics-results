library(magrittr)

metatable <- read.csv("../data/metatable_with_viral_status.csv",
                      stringsAsFactors=F)

probs_1se <- read.csv("../results/classifier/lasso_1se/lassoRandomForest_probs.csv",
                      stringsAsFactors=F)

probs_s10 <- read.csv("../results/classifier/lasso_s10/lassoRandomForest_probs.csv",
                      stringsAsFactors=F)

probs_s3 <- read.csv("../results/classifier/lasso_s3/lassoRandomForest_probs.csv",
                      stringsAsFactors=F)


get_sensitivity <- function(truth, pred) {
  mean((truth == pred)[truth])
}

get_specificity <- function(truth, pred) {
  mean((truth==pred)[!truth])
}

get_accuracy <- function(truth, pred) {
  mean(truth == pred)
}

get_ppv <- function(truth, pred) {
  mean((truth == pred)[pred])
}

get_npv <- function(truth, pred) {
  mean((truth == pred)[!pred])
}

get_subtable <- function(name, probs, threshold) {
  metatable %>%
    dplyr::select(CZB_ID, viral_status) %>%
    dplyr::inner_join(probs) %>%
    dplyr::mutate(truth=viral_status=="SC2",
                  pred=COVID19 >= threshold) %>%
    with(data.frame(
      Name=name,
      Threshold=threshold,
      Accuracy=get_accuracy(truth, pred),
      Sensitivity=get_sensitivity(truth, pred),
      Specificity=get_specificity(truth, pred),
      `Positive Predictive Value`=get_ppv(truth, pred),
      `Negative Predictive Value`=get_npv(truth, pred)
    ))
}

rbind(
  get_subtable("26-gene model", probs_1se, 0.5),
  get_subtable("10-gene model", probs_s10, 0.5),
  get_subtable("3-gene model", probs_s3, 0.5),
  get_subtable("26-gene model", probs_1se, 0.4),
  get_subtable("10-gene model", probs_s10, 0.4),
  get_subtable("3-gene model", probs_s3, 0.4)
) %>%
  format(digits=3) %>%
  write.csv("../figures/classifier/accuracy-table.csv", row.names=F)
