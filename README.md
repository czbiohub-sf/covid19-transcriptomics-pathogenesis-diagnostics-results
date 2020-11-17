# covid19-transcriptomics-pathogenesis-diagnostics-results
This repo houses data and code for the analyses in the paper ["Upper airway gene expression reveals suppressed immune responses to SARS-CoV-2 compared with other respiratory viruses"](https://doi.org/10.1038/s41467-020-19587-y) by Mick, Kamm, et al. published in *Nature Communications*.<br> [Note: the code associated with the [pre-print](https://doi.org/10.1101/2020.05.18.20105171) initially posted on medRxiv remains available for download as a release.]  

The analyses include:
1. Gene differential expression (DE) on metagenomic RNA-sequencing data from nasopharyngeal/oropharyngeal swabs of patients with COVID-19, other viral acute respiratory illnesses (ARIs) or non-viral ARIs.
2. Cell type deconvolution of the bulk RNA-seq in the three patient groups.
3. Diagnostic classifier based on patient gene expression to distinguish COVID-19 from other ARIs (viral or non-viral).

## Sample metadata table

Sample metadata is at `data/metatable.csv`. This table has the following fields:
1. `CZB_ID` - unique sample identifier.
2. `Sequencing_batch` - one of two library preparation and sequencing batches in which the samples were processed.
3. `Gender` - imputed based on chr Y gene expression.
4. `Age` - as reported by the patient (missing for some). Patient age was censored at 89 to comply with privacy requirements but this has negligible effects on the results. 
5. `SC2_PCR` - the result of the clinical PCR test for SARS-CoV-2. This field determines whether a patient belongs to the COVID-19 patient group.
6. `SC2_rpm` - reads-per-million mapping to SARS-CoV-2 in the metagenomic sequencing data. These values were generated using the pipeline here: https://github.com/czbiohub/sc2-msspe-bioinfo. The human subtracted FASTQ files are available for [download](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA633853).
7. `IDSeq_sample_name` - name of the matching IDSeq sample (see below).<br><br>
A second table `data/metatable_with_viral_status.csv` is generated by `code/viral_calls-2.R`, which determines the final assignment of samples with a negative SARS-CoV-2 PCR test either to the non-viral ARI or other viral ARI groups, depending on whether other pathogenic respiratory viruses were detected in the metagenomic sequencing (see below). This table includes an additional field:
8. `viral_status` - either SC2 (for SARS-CoV-2), other_virus or no_virus.<br>
Furthermore, the age field in this table imputes the age for patients for whom we lacked this information as the mean age of patients in their viral status group. 

## Viral calls

Viral calls are generated by the script `code/viral_calls-1.R`, which
creates `results/viral_calls/viruses_with_pval.csv`, a table
of viral read counts, along with p-values for being above background
levels according to a negative binomial model. The second script
`code/viral_calls-2.R` aggregates these results and combines them with
the metadata table to create `data/metatable_with_viral_status.csv`.

This script relies on code from the package
[idseqr](https://github.com/czbiohub/idseqr) which is currently in
alpha stage. The specific version of that package used for the
analysis is included as a git submodule at `code/idseqr`.

The script expects the metagenomic count data to reside in
`data/viral_calls/sample_taxon_reports`, but this data is not included in
the repo due to size limitations. To download the data, register for a
free account on [IDSeq](http://www.idseq.net), then download the "Sample Taxon
Reports" for the project
`covid19_transcriptomics_pathogenesis_diagnostics`. The Download will
prompt you to select a "Background" but it doesn't affect any of the
variables we use so you can select any background (e.g., "NID Human
CSF"). Unzip the reports and then move them into
`data/viral_calls/sample_taxon_reports`. The metadata table described above includes the IDSeq sample name associated with each sample CZB_ID.

## Gene counts table

Transcript quantification was performed with kallisto (v. 0.46.1, with bias correction) against an index consisting of transcripts of protein coding genes (ENSEMBL v. 99), cytosolic and mitochondrial ribosomal RNA and ERCC RNA standards.<br><br>
Aggregation to the gene level was done with tximport using countsFromAbundance="lengthScaledTPM". The mapping from transcript ID to gene ID is in `annotations/tx2gene.txt`. Genes were retained for analysis if they had at least 10 counts in at least 20% of the samples in the dataset. The counts for those gene are in `data/swab_gene_counts.csv`. Note, gene counts in the table were not normalized further. Genes in the table are identified by their ENSEMBL ID. The mapping to gene symbol and chromosome is in `annotations/gene2name.txt`.

## DE analysis

Pairwise DE analyses between the three patient groups using the model ~viral_status + age + gender, as well as regression of gene counts against viral abundance (rpM) for the differentially expressed genes, are performed in `code/DE_analysis.R`. Results are in `results/DE`. 

## Classifier analysis

- `classifier-1-generate-cv-folds.R` generates the
  cross-validation folds.
- `classifier-2-fit-models.R` fits the classifier models.
  Results are stored in `results/classifier`.
- `classifier-3-aggregate-results.py` and
  `classifier-4-sensitivity-specifity.R` aggregate results and put
  them in `figures/classifier`.

## Deconvolution analysis
We used the Human Lung Cell Atlas dataset [Travaglini et al. (bioRxiv 2019)](https://www.biorxiv.org/content/10.1101/742320v1) to generate the single cell signature files used as the basis for the deconvolution analysis in [CiberSortX](https://www.nature.com/articles/s41587-019-0114-2).

Cell deconvolution analysis is in `code/host_deconvolution.ipynb`,
with results in `results/deconvolution`.
