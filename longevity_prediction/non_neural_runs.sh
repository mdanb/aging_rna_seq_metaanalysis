python3 cv_non_neural.py --expression_path ../common_datastore/getmm_gene_expression_no_outliers.csv --mixsplit --train_XGB

python3 cv_non_neural.py --expression_path ../common_datastore/getmm_gene_expression_no_outliers.csv --train_XGB
python3 cv_non_neural.py --expression_path ../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv --train_XGB
python3 cv_non_neural.py --expression_path ../common_datastore/combat_seq_age_corrected_getmm_gene_expression_no_outliers.csv
python3 cv_non_neural.py --expression_path ../common_datastore/combat_seq_experiment_and_age_corrected_getmm_gene_expression_no_outliers.csv

python3 cv_non_neural.py --expression_path ../common_datastore/SAUCIE_getmm.csv

python3 cv_non_neural.py --aging_genes_only --expression_path ../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv
python3 cv_non_neural.py --expression_path ../common_datastore/combat_seq_experiment_corrected_day_1_and_older_getmm_gene_expression_no_outliers.csv
python3 cv_non_neural.py --expression_path ../common_datastore/combat_seq_experiment_corrected_L4_and_younger_getmm_gene_expression_no_outliers.csv
python3 cv_non_neural.py --expression_path ../common_datastore/combat_seq_getmm_GO_filtered_gene_expression_no_singles_and_outliers.csv
