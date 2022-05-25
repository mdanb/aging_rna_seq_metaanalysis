from joblib import load
import argparse
from sklearn.metrics import accuracy_score
from utils import load_data, load_train_test_datasets

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

dataset = load_data(config)
train_test_dataset_list = load_train_test_datasets(dataset)
X_test, labels_test, _ = train_test_dataset_list[1]

preds = best_clf.predict(X_test.values)
accuracy = accuracy_score(preds, labels_test)
with open("results/non_neural/best_model_test_set_performance.txt", "w") as f:
    f.write(f"{accuracy}\n")