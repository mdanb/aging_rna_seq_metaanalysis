from joblib import load
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('--expression_path', type=str,
                    default="../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv",
                    help='path to gene expression data '
                         '(default: ../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv)')
parser.add_argument('--label_path', type=str, default="../common_datastore/labels.csv",
                    help='path to labels (default: ../common_datastore/labels.csv)')
parser.add_argument('--age_path', type=str, default="../common_datastore/age.csv",
                    help='path to age data (default: ../common_datastore/age.csv)')
parser.add_argument('--experiments_path', type=str, default="../common_datastore/sra_to_bioproject.csv",
                    help='path to sra to bioproject mapping (default: ../common_datastore/sra_to_bioproject.csv)')
parser.add_argument('--aging_genes_only', action='store_true',
                    help="train the model using aging genes only (default: False)")


config = parser.parse_args()

aging_genes_only = config.aging_genes_only
data_filename = config.expression_path.split('/')[-1]
best_model = f"results/non_neural/aging_genes_only_{aging_genes_only}_{data_filename}_" \
             f"trained_best_non_neural_model.joblib"
best_clf = load(best_model)

longlived_sorted_features = np.argsort(best_clf[0].selected_model.coef_[0])[::-1]
genes = pd.read_csv("../common_datastore/getmm_combat_seq_no_outliers_and_singles_"
                      "gene_expression.csv", nrows=1).columns.tolist()[1:]
features = genes + ['age']
top_ten_long_idxs = longlived_sorted_features[:10]
top_long_features = [features[i] for i in top_ten_long_idxs]
bottom_ten_long_idxs = longlived_sorted_features[-10:]
bottom_long_features = [features[i] for i in bottom_ten_long_idxs]
print(top_long_features)
shortlived_sorted_features = np.argsort(best_clf[0].selected_model.coef_[2])[::-1]
top_ten_short_idxs = shortlived_sorted_features[:10]
top_short_features = [features[i] for i in top_ten_short_idxs]
bottom_ten_short_idxs = shortlived_sorted_features[-10:]
bottom_short_features = [features[i] for i in bottom_ten_short_idxs]
