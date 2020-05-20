library(magrittr)

rep_viruses_txt <- "../annotation/known_respiratory_pathogens.txt"
viral_calls_csv <- "../results/viral_calls/viruses_with_pval.csv"

metatable_csv <- "~/Box Sync/COVID-mNGS-Study/COVID-host-expression/metatable_age_uncensored_correct_idseq.csv"
out_table_csv <- "~/Box Sync/COVID-mNGS-Study/COVID-host-expression/metatable_age_uncensored_with_viral_status.csv"
## UNCOMMENT TO USE PUBLIC (AGE-CENSORED) DATA
#metatable_csv <- "../data/metatable.csv"
#out_table_csv <- "../data/metatable_with_viral_status.csv"

resp_viruses <- read.table(rep_viruses_txt,
                           header = FALSE, sep="\t", quote = "", col.names = "Pathogen")
resp_viruses %>%
  dplyr::filter(grepl("virus",Pathogen) & !grepl("mastadenovirus|betaherpesvirus",Pathogen,ignore.case = T)) %>%
  tibble::add_row(Pathogen = "Wuhan seafood market pneumonia virus") ->
  resp_viruses

metatable <- read.csv(metatable_csv)

viral_calls <- read.csv(viral_calls_csv, header = TRUE)

viral_calls %>%
  dplyr::mutate(sample_name=as.character(sample_name)) %>%
  dplyr::filter(sample_name %in% metatable$IDSeq_sample_name) %>%
  dplyr::filter(name %in% resp_viruses$Pathogen | common_name %in% resp_viruses$Pathogen,
         nt_count >= 10,
         nt_alignment_length >= 70,
         nr_count >= 1) %>%
  dplyr::group_by(sample_name) %>%
  dplyr::mutate("p_adj" = p.adjust(p_val, method = "holm")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(p_adj < 0.05) ->
  viral_calls

metatable %>%
  dplyr::mutate(
    viral_status=dplyr::case_when(
      SC2_PCR=="POS" ~ "SC2",
      IDSeq_sample_name %in% dplyr::filter(
        viral_calls, !startsWith(name, "Wuhan"))$sample_name ~ "other_virus",
      TRUE ~ "no_virus")) %>%
  write.csv(out_table_csv, row.names=F, quote=F)
