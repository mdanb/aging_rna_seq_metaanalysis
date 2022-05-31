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

python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --dropout 0.2
python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --dropout 0.5
python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --dropout 0.8

python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --weight_decay 0.1
python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --weight_decay 0.01
python3 cv_neural.py --mlp_hidden_dim 1024 --learning_rate 0.0001 --num_mlp_layers 3 --train_MLP --weight_decay 0.001
