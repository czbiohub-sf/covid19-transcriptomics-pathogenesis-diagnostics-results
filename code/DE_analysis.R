library(dplyr)
library(tibble)
library(magrittr)
library(limma)
library(robustbase)

metatable_csv <- "Box Sync/COVID-host-expression/metatable_age_uncensored_with_viral_status.csv"
## UNCOMMENT TO USE PUBLIC (AGE-CENSORED) DATA
#metatable_csv <- "../data/metatable_with_viral_status.csv"

counts_csv <- "../data/swab_gene_counts.csv"
gene_annot_tsv <- "../annotation/gene2name.txt"

# Read data
metatable_csv %>%
  read.csv() ->
  metatable

counts_csv %>%
  read.csv(row.names=1) %>% 
  as.matrix() ->
  counts

gene_annot_tsv %>%
  read.delim(header = F, col.names = c("gene_ID","gene_symbol", "chr"), row.names=1) ->
  gene_annot

metatable$viral_status <- factor(metatable$viral_status, levels = c("SC2","no_virus","other_virus"))

# Set up design matrix
design_matrix <- model.matrix(~0 + viral_status + gender + age + sequencing_batch, data = metatable)
colnames(design_matrix) <- c("SC2", "no_virus", "other_virus", "gender_M", "age", "SEQ003")

# limma-voom
vwts <- voom(counts, 
             design = design_matrix,
             normalize.method = "quantile") 
vfit <- lmFit(vwts)

# Set up contrasts
contr <- makeContrasts(SC2-no_virus, 
                       SC2-other_virus,
                       other_virus-no_virus,
                       levels = vfit$design)
contr_fit <- contrasts.fit(vfit, contr)
contr_fit <- eBayes(contr_fit)

# Extract contrasts
sc2_vs_no <- topTable(contr_fit, coef = 1, sort.by = "none", number = nrow(counts))
sc2_vs_other <- topTable(contr_fit, coef = 2, sort.by = "none", number = nrow(counts))
other_vs_no <- topTable(contr_fit, coef = 3, sort.by = "none", number = nrow(counts))

# Combine DE results into one dataframe
combined_de <- data.frame("gene_symbol" = gene_annot[rownames(sc2_vs_no), "gene_symbol"],
                          "SC2_no_virus_padj" = sc2_vs_no$adj.P.Val,
                          "SC2_no_virus_logFC" = sc2_vs_no$logFC,
                          "SC2_other_virus_padj" = sc2_vs_other$adj.P.Val,
                          "SC2_other_virus_logFC" = sc2_vs_other$logFC,
                          "other_virus_no_virus_padj" = other_vs_no$adj.P.Val,
                          "other_virus_no_virus_logFC" = other_vs_no$logFC,
                          row.names = rownames(sc2_vs_no))

# Pick genes to regress on SC2 rpm
p_adj_thresh <- 1e-3
genes_to_reg <- union(rownames(sc2_vs_no[sc2_vs_no$adj.P.Val < p_adj_thresh,]),
                      rownames(sc2_vs_other[sc2_vs_other$adj.P.Val < p_adj_thresh,]))

# Limit to samples that are positive for SARS-CoV-2 by PCR and have at least 1 rpm by mNGS
samples_to_reg <- metatable %>% filter(viral_status=="SC2", SC2_rpm >= 1) %>% pull(CZB_ID)

# Perform robust regression of the limma-generated quantile normalized counts (log2 scale) on log10(SC2_rpm + 0.1); +0.1 required to avoid log(0)
rpm_gene_robs <- as.data.frame(t(sapply(genes_to_reg, function(g) {
  reg_data <- metatable %>% mutate("log_norm_counts" = vwts$E[g, ], "log_rpm" = log10(SC2_rpm + 0.1)) %>% filter(CZB_ID %in% samples_to_reg) 
  rob_reg <- summary(lmrob(log_norm_counts ~ gender + age + sequencing_batch + log_rpm, 
                           data = reg_data,
                           setting = "KS2014"))
  c("reg_intercept" = rob_reg$coefficients[1,1], 
    "reg_slope" = rob_reg$coefficients[5,1], 
    "reg_p_val" = rob_reg$coefficients[5,4], 
    "reg_adj_R2" = rob_reg$adj.r.squared)
})))

# Adjust regression p-values
rpm_gene_robs[,"reg_p_adj"] <- p.adjust(rpm_gene_robs$reg_p_val, method="BH")

# Combine regression results into DE table
combined_de %>% 
  rownames_to_column("gene_id") %>% 
  left_join(x = ., y = rpm_gene_robs %>% rownames_to_column("gene_id"), by = "gene_id") %>% 
  column_to_rownames("gene_id") %>% 
  select(-reg_p_val) ->
  combined_de

# Write DE+regression table
write.csv(combined_de, "../results/DE/3way_diff_expr_and_regression.csv")

