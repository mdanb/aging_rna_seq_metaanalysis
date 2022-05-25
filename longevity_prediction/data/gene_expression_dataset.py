from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import data.utils

class GeneDataset(Dataset):
    def __init__(self):
        self.load_data()

    def load_data(self):
        raise NotImplementedError()

    def __getitem__(self, idx):
        raise NotImplementedError()


class DatasetFromCSV(GeneDataset):
    def __init__(self, expr_path, label_path, age_path, experiments_path, aging_genes_only):
        self.expr_path = expr_path
        self.label_path = label_path
        self.age_path = age_path
        self.experiments_path = experiments_path
        self.aging_genes_only = aging_genes_only
        super(DatasetFromCSV, self).__init__()

    def load_data(self):
        sep = data.utils.get_file_separator(self.expr_path)
        self.df = pd.read_csv(self.expr_path, sep=sep, index_col=0)
        if (self.aging_genes_only):
            sep = data.utils.get_file_separator('../../common_datastore/genage_gene_name_to_wb_id.txt')
            aging_genes = pd.read_csv("../common_datastore/genage_gene_name_to_wb_id.txt", sep=sep).iloc[:, 1]
            aging_genes = [gene for gene in aging_genes if gene in self.df.columns]
            self.df = self.df[aging_genes]
        sep = data.utils.get_file_separator(self.label_path)
        self.labels = pd.read_csv(self.label_path, sep=sep, index_col=0)
        self.labels = self.df.merge(self.labels, left_index=True, right_index=True).iloc[:, -1]
        self.labels = pd.DataFrame(self.labels)
        self.labels.columns = ['longevity']

        sep = data.utils.get_file_separator(self.age_path)
        self.age = pd.read_csv(self.age_path, sep=sep, index_col=0)
        self.age = self.df.merge(self.age, left_index=True, right_index=True).iloc[:, -1]
        self.age = pd.DataFrame(self.age)
        self.age.columns = ['age']

        sep = data.utils.get_file_separator(self.experiments_path)
        self.experiments = pd.read_csv(self.experiments_path, sep=sep, index_col=0)['Bioproject']
        self.experiments = self.df.merge(self.experiments, left_index=True, right_index=True).iloc[:, -1]

    def __getitem__(self, idx):
        sample = self.df.iloc[idx,:].values
        sample = np.expand_dims(sample, axis=-1)
        label = self.labels.values[idx]
        age = self.age.values[idx]
        sample = {'x': sample, 'y': label, 'age': age}
        return sample

    def __len__(self):
        pass

