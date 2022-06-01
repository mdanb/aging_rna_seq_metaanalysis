Here, you will find the scripts that were used to Combat-Seq/GeTMM normalize the data, as well as those used to create the PCA plots and similarity matrices.

To Combat-Seq/GeTMM the data, run `Rscript getmm_and_combat_seq.R`

To create PCA plots and similarity matrices, run `python eda.py` (the contents of `eda.ipynb` are the same as those in `eda.py`, and were used to run the code in a Jupyter notebook). 

Note: for anything GO related, you will first need to run `python ../longevity_prediction/find_filtered_genes.py`
