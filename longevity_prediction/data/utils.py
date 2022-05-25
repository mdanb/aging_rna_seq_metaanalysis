import pickle
from gtfparse import read_gtf
import os


def get_file_separator(filename):
    separators = {'.tsv' : '\t', '.txt': '\t', '.csv': ','}
    sep = separators[os.path.splitext(filename.replace('.gz', ''))[-1]]
    return sep

def ensembl_protein_to_gene_map():
    savefile = "./data/datastore/ensembl_protein_gene_map.pkl"

    if os.path.isfile(savefile):
        f = open(savefile, 'rb')
        df = pickle.load(f)
        f.close()
    else:
        df = read_gtf("./data/datastore/Caenorhabditis_elegans.WBcel235.100.gtf")
        df = df[df['protein_id'] != ''][['gene_id', 'protein_id']].drop_duplicates()
        df.to_pickle(savefile)

    # return ENSP to ENSG dictionary
    return {row['protein_id']: row['gene_id'] for index, row in df.iterrows()}