library(magrittr)

seed = 210274353
metatable_csv = "../data/metatable_with_mngs_viral_calls.csv"
folds_out_csv = "../results/cv-folds.csv"
log_out_csv = "../data/cv-folds.log"

## Set seed

set.seed(seed)

## Read data

metatable_csv %>%
  readr::read_csv() %>%
  as.data.frame() ->
  metatable

## Create folds. Balance class sizes due to relatively small size of
## OtherVirus.

covid_idxs <- metatable$test_result_final=="POS"
other_idxs <- !covid_idxs & metatable$Other_virus_by_mNGS
no_virus_idxs <- !covid_idxs & !other_idxs

folds <- rep(0, length.out=nrow(metatable))
folds[covid_idxs] <- sample(rep(1:5, length.out=sum(covid_idxs)))
folds[other_idxs] <- sample(rep(1:5, length.out=sum(other_idxs)))
folds[no_virus_idxs] <- sample(rep(1:5, length.out=sum(no_virus_idxs)))

stopifnot(folds != 0)

## Save folds

metatable %>%
  dplyr::select(CZB_ID) %>%
  dplyr::mutate(fold=folds) %>%
  write.csv(folds_out_csv, row.names=F)

## Save version info and seed

sink(file=log_out_csv)
cat(paste("Seed: " , seed, "\n\n", sep=""))
sessionInfo()
sink(file=NULL)
