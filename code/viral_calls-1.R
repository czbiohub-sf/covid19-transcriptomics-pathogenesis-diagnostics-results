library(magrittr)
library(ggplot2)
library(idseqr)

### paths

idseq_sample_overviews_path = "../data/viral_calls/sample_overviews.csv"
ngs_samples_path = "../data/viral_calls/ngs_samples.csv"
known_pathogens_path = "../annotation/known_respiratory_pathogens.txt"
reports_path = "../data/viral_calls/sample_taxon_reports"
results_xy_path = "../results/viral_calls/rpm_vs_input.svg"
results_heatmap_path = "../results/viral_calls/heatmap.svg"
results_pval_path = "../results/viral_calls/viruses_with_pval.csv"

### Read data

read.csv(idseq_sample_overviews_path, stringsAsFactors=F) %>%
  dplyr::mutate(
    control=dplyr::case_when(
      water_control=="Yes" ~ TRUE,
      grepl("hela", sample_name, ignore.case=T) ~ TRUE,
      TRUE ~ FALSE)) %>%
  dplyr::inner_join(read.csv(ngs_samples_path, stringsAsFactors=F)) ->
  sample_overviews

known_resp <- readLines(known_pathogens_path)

read_reports(reports_path, tax_level=2) %>%
  dplyr::filter(sample_name %in% sample_overviews$sample_name) ->
  genus_reports

read_reports(reports_path) %>%
  dplyr::filter(sample_name %in% sample_overviews$sample_name) %>%
  dplyr::mutate(genus_tax_id=dplyr::if_else(genus_tax_id < 0,
                                            tax_id, genus_tax_id)) %>%
  dplyr::left_join(
    genus_reports %>%
      dplyr::select(sample_name, tax_id, nt_count) %>%
      dplyr::rename(genus_tax_id=tax_id, genus_nt_count=nt_count)
  ) %>%
  dplyr::mutate(genus_nt_count=dplyr::if_else(is.na(genus_nt_count),
                                              nt_count, genus_nt_count)) ->
  reports

stopifnot(length(unique(reports$sample_name)) == nrow(sample_overviews))
stopifnot(length(unique(genus_reports$sample_name)) == nrow(sample_overviews))

### Background correction

controls <- with(sample_overviews, sample_name[control])

ercc_norm <- sample_overviews$total_ercc_reads
names(ercc_norm) <- sample_overviews$sample_name

batches <- sample_overviews$sequencing_batch
names(batches) <- sample_overviews$sample_name

reports_filterbg <- filter_background(reports, controls,
                                      batches=batches,
                                      normalization=ercc_norm)

### Plots

svg(results_xy_path, width=12, height=8)
reports_filterbg %>%
  #filter_top_taxa(top_tax_per_sample=1) %>%
  dplyr::filter(category=="viruses") %>%
  dplyr::group_by(sample_name) %>%
  dplyr::mutate(p_adj=p.adjust(p_val, method="holm")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(name %in% c("Human mastadenovirus C",
                            "Human orthopneumovirus",
                            "Wuhan seafood market pneumonia virus")) %>%
  dplyr::inner_join(sample_overviews) %>%
  ggplot(aes(x=(total_reads-total_ercc_reads)/total_ercc_reads, y=nt_rpm,
             pch=p_adj <= .05, color=control)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  scale_shape_manual(values=c(21, 8)) +
  facet_wrap(~name+sequencing_batch)
dev.off()

svg(results_heatmap_path, width=20, height=8)
reports_filterbg %>%
  dplyr::filter(category=="viruses") %>%
  dplyr::filter(name %in% known_resp) %>%
  dplyr::filter(nt_count / genus_nt_count >= .1, nt_alignment_length >= 70) %>%
  dplyr::group_by(sample_name) %>%
  dplyr::mutate(p_adj=p.adjust(p_val, method="holm")) %>%
  dplyr::ungroup() %>%
  as.data.frame() %>%
  filter_top_taxa(1, counts_column="nt_count") %>%
  dplyr::inner_join(sample_overviews) %>%
  {
    df <- .
    ggplot(df, aes(x=sample_name, y=name, fill=nt_count)) +
      geom_tile() +
      geom_point(data=dplyr::filter(df, p_adj < .05, nt_count >= 5)) +
      facet_grid(category ~ sequencing_batch + control,
                 scales="free", space="free", labeller="label_both") +
      scale_fill_gradientn(colours=c("turquoise", "violet", "orange"), trans="log10") +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1))
  }
dev.off()

reports_filterbg %>%
  dplyr::filter(category=="viruses") %>%
  write.csv(results_pval_path, row.names=F)
