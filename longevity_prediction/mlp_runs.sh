python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.001 --num_mlp_layers 5 --train_MLP
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.01 --num_mlp_layers 5 --train_MLP 
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 5 --train_MLP

python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 4 --train_MLP
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 6 --train_MLP
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 2 --train_MLP

python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP
python3 cv_neural.py --mlp_hidden_dim 256 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP
python3 cv_neural.py --mlp_hidden_dim 2048 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP

python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --dropout 0.2
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --dropout 0.5
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --dropout 0.8

python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --weight_decay 0.1
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --weight_decay 0.01
python3 cv_neural.py --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --weight_decay 0.001

#python3 cv_neural.py --expression_path ../common_datastore/combat_seq_age_corrected_day_1_and_older_getmm_gene_expression_no_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP
#python3 cv_neural.py --expression_path ../common_datastore/combat_seq_age_corrected_L4_and_younger_getmm_gene_expression_no_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP
#python3 cv_neural.py --expression_path ../common_datastore/combat_seq_getmm_GO_filtered_gene_expression_no_singles_and_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP
#python3 cv_neural.py --expression_path ../common_datastore/combat_seq_age_corrected_day_1_and_older_getmm_gene_expression_no_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP
#python3 cv_neural.py --expression_path ../common_datastore/combat_seq_experiment_and_age_corrected_getmm_gene_expression_no_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP
#python3 cv_neural.py --expression_path ../common_datastore/double_combat_seq_experiment_corrected_getmm_gene_expression_no_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP
#python3 cv_neural.py --expression_path ../common_datastore/getmm_gene_expression_no_outliers.csv --mlp_hidden_dim 512 --learning_rate 0.0001 --num_mlp_layers 3 --weight_decay 0.001 --train_MLP

