import pandas as pd
from sklearn.decomposition import PCA
from utils import load_data, load_train_test_datasets
import argparse
from split import StratifiedGroupKFold
import os

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

def reduce_dim_to_minimum_variance(X_train, labels_train,
                                   experiments_train,
                                   minimum_variances):
    for minimum_variance in minimum_variances:
        cv = StratifiedGroupKFold(10)
        reduced_dim_X_train_val = []
        for train_idxs, val_idxs in cv.split(X_train,
                                             labels_train.to_numpy().reshape(-1),
                                             experiments_train):
            X_train_cv = X_train.iloc[train_idxs, :]
            X_val_cv = X_train.iloc[val_idxs, :]
            pca = PCA(n_components=minimum_variance)
            fit_pca = pca.fit(X_train_cv)
            X_val_cv = fit_pca.transform(X_val_cv)

config = parser.parse_args()
dataset = load_data(config)
train_test_dataset_list = load_train_test_datasets(dataset)
X_train, labels_train, experiments_train = train_test_dataset_list[0]
reduce_dim_to_minimum_variance(X_train, labels_train,
                               experiments_train,
                               [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99])
