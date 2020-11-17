library(magrittr)

rep_viruses_txt <- "../annotation/known_respiratory_pathogens.txt"
viral_calls_csv <- "../results/viral_calls/viruses_with_pval.csv"

metatable_csv <- "../data/metatable.csv"
out_table_csv <- "../data/metatable_with_viral_status.csv"

resp_viruses <- read.table(rep_viruses_txt,
                           header = FALSE, sep="\t", quote = "", col.names = "Pathogen")
resp_viruses %>%
  dplyr::filter(grepl("virus",Pathogen) & !grepl("mastadenovirus|betaherpesvirus",Pathogen,ignore.case = T)) %>%
  tibble::add_row(Pathogen = "Wuhan seafood market pneumonia virus") ->
  resp_viruses

metatable <- read.csv(metatable_csv, stringsAsFactors = F)

viral_calls <- read.csv(viral_calls_csv, header = TRUE, stringsAsFactors = F)

viral_calls %>%
  dplyr::mutate(sample_name=as.character(sample_name)) %>%
  dplyr::filter(sample_name %in% metatable$idseq_sample_name) %>%
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
      idseq_sample_name %in% dplyr::filter(
        viral_calls, !startsWith(name, "Wuhan"))$sample_name ~ "other_virus",
      TRUE ~ "no_virus")) ->
  metatable

# Calculate mean age in each viral status group
metatable %>% 
  dplyr::group_by(viral_status) %>% 
  dplyr::summarize("mean_age" = mean(age, na.rm = T)) %>% 
  tibble::column_to_rownames("viral_status") ->
  mean_age_by_viral_status

# Set the age of patients with missing age to the mean age in their viral status group
metatable %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate("age" = ifelse(is.na(age),
                        mean_age_by_viral_status[viral_status, "mean_age"],
                        age)) %>% 
  dplyr::ungroup() ->
  metatable

metatable %>%
  write.csv(out_table_csv, row.names=F, quote=F)