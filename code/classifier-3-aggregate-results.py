import os
import csv
import numpy as np
import pandas as pd

def read_auc_file(auc_file):
    auc_arr = []
    with open(auc_file) as f:
        for line in f:
            auc_arr.append(float(line.strip()))
    auc_arr = np.array(auc_arr)
    return "{:.3f} ({:.3f}-{:.3f})".format(auc_arr.mean(), auc_arr.min(), auc_arr.max())

def read_nonzero_coefs_table(path):
    ret = pd.read_csv(path)
    corrected_gene_names = []
    for _, row in ret.iterrows():
        if pd.isnull(row["gene_name"]):
            corrected_gene_names.append(row["variable"])
        else:
            corrected_gene_names.append(row["gene_name"])
    ret["gene_name"] = corrected_gene_names
    ret.drop("variable", axis=1, inplace=True)
    return ret

def read_results_folder(results_folder):
    return {
        "lasso": {
            "All": read_auc_file(
                os.path.join(results_folder, "lasso_cv_auc.txt")),
            "COVID vs No virus": read_auc_file(
                os.path.join(results_folder, "lasso_excludeOther_cv_auc.txt")),
            "Covid vs Other virus": read_auc_file(
                os.path.join(results_folder, "lasso_viralOnly_cv_auc.txt")),
        },
        "lasso + randomForest": {
            "All": read_auc_file(
                os.path.join(results_folder, "lassoRandomForest_cv_auc.txt")),
            "COVID vs No virus": read_auc_file(
                os.path.join(results_folder, "lassoRandomForest_excludeOther_cv_auc.txt")),
            "Covid vs Other virus": read_auc_file(
                os.path.join(results_folder, "lassoRandomForest_viralOnly_cv_auc.txt")),
        },
        "lasso_nonzero_coefs": read_nonzero_coefs_table(
            os.path.join(results_folder, "lasso_nonzero_coefs.csv"))
    }

results_dir = "../results/classifier"
out_dir = "../figures/classifier"

with open(os.path.join(out_dir, "auc_table.csv"), "w") as f:
    writer = csv.DictWriter(
        f, fieldnames=["Model", "All", "COVID vs No virus", "Covid vs Other virus"])

    writer.writeheader()

    lasso_1se = read_results_folder(os.path.join(results_dir, "lasso_1se"))

    row = lasso_1se["lasso + randomForest"]
    row["Model"] = "Full 26-gene model (lasso + randomForest)"
    writer.writerow(row)

    row = lasso_1se["lasso"]
    row["Model"] = "Full 26-gene model (lasso only)"
    writer.writerow(row)

    lasso_1se["lasso_nonzero_coefs"].to_csv(
        os.path.join(out_dir, "lasso_nonzero_coefs_1se.csv"),
        float_format="%.3f", index=False)

    lasso_s10 = read_results_folder(os.path.join(results_dir, "lasso_s10"))
    row = lasso_s10["lasso + randomForest"]
    row["Model"] = "10-gene model (lasso + randomForest)"
    writer.writerow(row)

    row = lasso_s10["lasso"]
    row["Model"] = "10-gene model (lasso only)"
    writer.writerow(row)

    lasso_s10["lasso_nonzero_coefs"].to_csv(
        os.path.join(out_dir, "lasso_nonzero_coefs_s10.csv"),
        float_format="%.3f", index=False)

    lasso_s5 = read_results_folder(os.path.join(results_dir, "lasso_s5"))
    row = lasso_s5["lasso + randomForest"]
    row["Model"] = "5-gene model (lasso + randomForest)"
    writer.writerow(row)

    row = lasso_s5["lasso"]
    row["Model"] = "5-gene model (lasso only)"
    writer.writerow(row)

    lasso_s5["lasso_nonzero_coefs"].to_csv(
        os.path.join(out_dir, "lasso_nonzero_coefs_s5.csv"),
        float_format="%.3f", index=False)

    lasso_s3 = read_results_folder(os.path.join(results_dir, "lasso_s3"))
    row = lasso_s3["lasso + randomForest"]
    row["Model"] = "3-gene model (lasso + randomForest)"
    writer.writerow(row)

    row = lasso_s3["lasso"]
    row["Model"] = "3-gene model (lasso only)"
    writer.writerow(row)

    lasso_s3["lasso_nonzero_coefs"].to_csv(
        os.path.join(out_dir, "lasso_nonzero_coefs_s3.csv"),
        float_format="%.3f", index=False)

    withAgeGender_1se = read_results_folder(os.path.join(results_dir, "withAgeGender_lasso_1se"))
    row = withAgeGender_1se["lasso + randomForest"]
    row["Model"] = "Model with Age and Gender covariates\n(20 genes, lasso + randomForest)"
    writer.writerow(row)

    withAgeGender_s10 = read_results_folder(os.path.join(results_dir, "withAgeGender_lasso_s10"))
    row = withAgeGender_s10["lasso + randomForest"]
    row["Model"] = "Model with Age and Gender covariates\n(10 genes, lasso + randomForest)"
    writer.writerow(row)

    withAgeGender_s5 = read_results_folder(os.path.join(results_dir, "withAgeGender_lasso_s5"))
    row = withAgeGender_s5["lasso + randomForest"]
    row["Model"] = "Model with Age and Gender covariates\n(5 genes, lasso + randomForest)"
    writer.writerow(row)

    withAgeGender_s3 = read_results_folder(os.path.join(results_dir, "withAgeGender_lasso_s3"))
    row = withAgeGender_s3["lasso + randomForest"]
    row["Model"] = "Model with Age and Gender covariates\n(3 genes, lasso + randomForest)"
    writer.writerow(row)
