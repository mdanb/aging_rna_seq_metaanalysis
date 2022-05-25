import pandas as pd
import os

if not (os.path.exists("../common_datastore/SAUCIE_getmm.csv")):
    df = pd.DataFrame()
    sra_to_bioproject = pd.read_csv("../common_datastore/sra_to_bioproject.csv")

    all_sras = []
    for bioproject in sra_to_bioproject.iloc[:, 1].unique():
        try:
            gene_expression = pd.read_csv(f"scripts/outputs/batch_corrected/gene_expression_{bioproject}.csv")
            sras = sra_to_bioproject[sra_to_bioproject['Bioproject'] == bioproject]['Sample']
            all_sras.extend(sras)
            df = pd.concat([df, gene_expression], axis=0)
        except:
            continue # outlier not considered
    df.index = all_sras
    df.columns = pd.read_csv('../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv',
                             nrows=1).columns[1:]
    df.to_csv("../common_datastore/SAUCIE_getmm.csv")
