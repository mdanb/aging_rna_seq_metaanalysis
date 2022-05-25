python3 cv_neural.py --learning_rate 0.001 --train_GNN
python3 cv_neural.py --learning_rate 0.01 --train_GNN
python3 cv_neural.py --learning_rate 0.0001 --train_GNN

python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 3 --train_GNN
python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 2 --train_GNN      
python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 1 --train_GNN

python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 3 --concat_input_graph --train_GNN      

python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 3 --concat_input_graph --k 15000 --batch_size 512 --train_GNN
python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 3 --concat_input_graph --k 10000 --batch_size 512 --train_GNN

python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 3 --num_backbone_layers 2 --batch_size 512 --concat_input_graph --train_GNN
python3 cv_neural.py --learning_rate 0.0001 --num_mlp_layers 3 --num_backbone_layers 3 --batch_size 512 --concat_input_graph --train_GNN
