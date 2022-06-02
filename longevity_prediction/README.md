To obtain non-neural cross validation results, run `sh non_neural_runs`. To obtain neural cross validation results, run `sh mlp_runs.sh` for the MLP and `sh gnn_runs.sh` for the GNN. To obtain test results for the best non-neural and neural models, run `python test_non_neural.py` and `python test_neural.py` respectively. 

All results will be saved in the `results` folder. 

To get the top predictor genes according to the best non-neural model, run `python interpret_best_model.py`
