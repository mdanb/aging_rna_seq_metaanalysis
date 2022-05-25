from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from split import StratifiedGroupKFold
from xgboost import XGBClassifier
from pipelinehelper import PipelineHelper
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV, StratifiedShuffleSplit
import os
import numpy as np
from joblib import dump
import argparse
from utils import load_data, load_train_test_datasets

parser = argparse.ArgumentParser()

# Gene Expression settings
# Options:
# getmm_gene_expression_no_outliers.csv for getmm, no outliers
# getmm_combat_seq_no_outliers_and_singles_gene_expression.csv for getmm, combat-seq, no singles, no outliers
# combat_seq_age_corrected_getmm_gene_expression_no_outliers.csv
# combat_seq_age_corrected_day_1_and_older_getmm_gene_expression_no_outliers.csv
# combat_seq_age_corrected_L4_and_younger_getmm_gene_expression_no_outliers.csv
# combat_seq_getmm_GO_filtered_gene_expression_no_singles_and_outliers.csv
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
# Training / Testing settings
parser.add_argument('--procs', type=int, default=4, help="number of processors for "
                                                         "performing GridSearch and training XGBoost")
parser.add_argument('--seed', type=int, default=42, help="Seed")
parser.add_argument('--mixsplit', action='store_true', help="perform a mixsplit as described "
                                                            "in paper")
parser.add_argument('--train_XGB', action='store_true')

os.makedirs('results', exist_ok=True)

def save_grid_search_prediction_statistics(cv, best_model, X, y, data_filename, experiments, stratify_by_experiment,
                                           aging_genes_only):
    for idx, (_, test_index) in enumerate(cv.split(X, y.values.reshape(-1), experiments)):
        if (stratify_by_experiment):
            preds_savefile = f"results/non_neural/group_stratified_best_model_predictions_fold_{idx}_" \
                             f"aging_genes_only_{aging_genes_only}_{data_filename}"
            pred_props_savefile = f"results/non_neural/group_stratified_best_model_predicted_proportions_fold_{idx}_" \
                                  f"aging_genes_only_{aging_genes_only}_{data_filename}"
        else:
            preds_savefile = f"results/non_neural/best_model_predictions_fold_{idx}_aging_genes_only_" \
                             f"{aging_genes_only}_{data_filename}"
            pred_props_savefile = f"results/non_neural/best_model_predicted_proportions_fold_{idx}_" \
                                  f"aging_genes_only_{aging_genes_only}_{data_filename}"

        predictions = pd.DataFrame(best_model.predict(X.iloc[test_index, :].values), index = X.index[test_index],
                                   columns=['Predictions'])
        predictions.to_csv(preds_savefile)
        predicted_proportions = pd.DataFrame(predictions.value_counts() / sum(predictions.value_counts()),
                                             columns= ['Proportions'])
        predicted_proportions.rename(index={0: "long-lived", 1: "normal-lived", 2: "short-lived"}, inplace=True)
        predicted_proportions.to_csv(pred_props_savefile)
        fold_data = pd.DataFrame(columns=['Sample', 'Experiment', 'Label'])
        fold_data['Sample'] = X.index[test_index]
        fold_data['Experiment'] = experiments.values[test_index]
        fold_data['Label'] = y.values[test_index]
        fold_data.to_csv(f"results/non_neural/fold_{idx}_data.csv", index=False)

def grid_search_and_save_results(X, y, pipe, params, stratify_by_experiment, data_filename, experiments,
                                 stats_savefile, aging_genes_only):
    if (stratify_by_experiment):
        print("Stratify by experiment and label")
        cv = StratifiedGroupKFold(10)
        grid = GridSearchCV(pipe, params, scoring='accuracy', n_jobs=config.procs, cv=cv, verbose=100,
                            error_score="raise")
        results = grid.fit(X=X.values, y=y.values.reshape(-1), groups=experiments)
    else:
        print("Stratify by label only")
        cv = StratifiedShuffleSplit(n_splits=10, test_size=0.1, random_state=config.seed)
        grid = GridSearchCV(pipe, params, scoring='accuracy', n_jobs=config.procs, cv=cv, verbose=100)
        results = grid.fit(X=X.values, y=y.values.reshape(-1))

    dump(results.best_estimator_, f"results/non_neural/aging_genes_only_{aging_genes_only}_{data_filename}_"
                                  f"trained_best_non_neural_model.joblib")
    save_grid_search_prediction_statistics(cv, results.best_estimator_, X, y, data_filename, experiments,
                                           stratify_by_experiment, aging_genes_only)
    statistics = pd.DataFrame(results.cv_results_)
    statistics.to_csv(stats_savefile, index=False)

def compute_per_fold_stats(stats_savefile):
    stats_df = pd.read_csv(stats_savefile)
    empty = np.empty([1, stats_df.shape[1]])
    empty[:] = np.NaN
    means_to_append = pd.DataFrame(empty)
    means = stats_df.iloc[:, 6:16].mean(axis=0)
    means_to_append.iloc[0, 6:16] = means
    means_to_append.columns = stats_df.columns
    stats_df = pd.concat([stats_df, means_to_append], axis=0)
    stds_to_append = pd.DataFrame(empty)
    stds = stats_df.iloc[:, 6:16].std(axis=0)
    stds_to_append.iloc[0, 6:16] = stds
    stds_to_append.columns = stats_df.columns
    stats_df = pd.concat([stats_df, stds_to_append], axis=0)
    stats_df.to_csv(stats_savefile, index=False)

config = parser.parse_args()

if (config.train_XGB):
    pipe = Pipeline([
        ('classifier', PipelineHelper([
            ('xgb', XGBClassifier(gpu_id=0)),
            ('rf', RandomForestClassifier()),
            ('svm', SVC(kernel='linear')),
            ('lr', LogisticRegression()),
        ])),
    ])
    params = {
        'classifier__selected_model': pipe.named_steps['classifier'].generate({
            'xgb__max_depth':[3, 6, 9],
            'rf__max_depth':[3, 6, 9],
            'svm__C': [0.1, 1.0, 10, 1000],
            'lr__C': [0.1, 1.0, 10, 1000],
        })
    }
else:
    pipe = Pipeline([
    ('classifier', PipelineHelper([
        ('rf', RandomForestClassifier()),
        ('svm', SVC(kernel='linear')),
        ('lr', LogisticRegression()),
    ])),
    ])
    params = {
        'classifier__selected_model': pipe.named_steps['classifier'].generate({
            'rf__max_depth':[3, 6, 9],
            'svm__C': [0.1, 1.0, 10, 1000],
            'lr__C': [0.1, 1.0, 10, 1000],
        })
    }
os.makedirs('results/non_neural', exist_ok=True)
data_filename = config.expression_path.split('/')[-1]
if (config.mixsplit):
    stats_savefile = f"results/non_neural/non_neural_statistics_aging_genes_only_{config.aging_genes_only}_" \
                     f"{data_filename}"
    if (not os.path.exists(stats_savefile)):
        dataset = load_data(config)
        train_test_dataset_list = load_train_test_datasets(dataset)
        X_train, labels_train, experiments_train = train_test_dataset_list[0]
        grid_search_and_save_results(pd.DataFrame(X_train), labels_train, pipe, params,
                                        False, data_filename,
                                        experiments_train, stats_savefile, config.aging_genes_only,
                                        )
        # Compute per fold stats
        compute_per_fold_stats(stats_savefile)
else:
    stats_savefile = f"results/non_neural/non_neural_statistics_group_stratified_aging_genes_only_" \
                     f"{config.aging_genes_only}_{data_filename}"
    if (not os.path.exists(stats_savefile)):
        dataset = load_data(config)
        train_test_dataset_list = load_train_test_datasets(dataset)
        X_train, labels_train, experiments_train = train_test_dataset_list[0]
        grid_search_and_save_results(pd.DataFrame(X_train), labels_train, pipe, params,
                                        True, data_filename,
                                        experiments_train, stats_savefile, config.aging_genes_only,
                                        )
        # Compute per fold stats
        compute_per_fold_stats(stats_savefile)