from torch_geometric.data import Data
from torch_geometric.data import Dataset
import torch
import os.path as osp

class GeneExpressionGraph(Dataset):
    def __init__(self, gene_interaction_graph, dataset, config, labels, age,
                 transform=None, pre_transform=None, pre_filter=None):
        self.config = config
        self.dataset = dataset
        self.gene_interaction_graph = gene_interaction_graph
        self.file_extension = config.expression_path.split('/')[-1]
        self.labels = labels
        self.age = age
        if (osp.exists(f"./data/datastore/gene_exp_graph_pre_exp_values_merge.pt")):
            self.gene_exp_graph_pre_exp_values_merge = torch.load(
                f"./data/datastore/gene_exp_graph_pre_exp_values_merge.pt")

        super(GeneExpressionGraph, self).__init__("./data/datastore", transform,
                                                  pre_transform, pre_filter)

    def process(self):
        print("Creating gene expression graphs")
        print("This may take a while the first time...")

        if (not osp.exists(f"./data/datastore/gene_exp_graph_pre_exp_values_merge.pt")):
            dataset_genes = self.dataset.columns.tolist() # Gene index in the list = gene id for the
                                                                # purposes of mapping
            dataset_gene_edges = []
            dataset_gene_edge_weights = []

            for idx, edge in enumerate(self.gene_interaction_graph.graph_edgelist):
                if edge[0] in dataset_genes and edge[1] in dataset_genes:
                    dataset_gene_edges.append([dataset_genes.index(edge[0]),
                                              dataset_genes.index(edge[1])])
                    dataset_gene_edge_weights.append(self.gene_interaction_graph.edge_weights[idx])

            edge_attr = dataset_gene_edge_weights

            gene_exp_graph_pre_exp_values_merge = Data(
                    edge_index=torch.tensor(dataset_gene_edges, dtype=torch.long).t().contiguous(),
                    edge_attr=edge_attr)

            torch.save(gene_exp_graph_pre_exp_values_merge,
                       f"./data/datastore/gene_exp_graph_pre_exp_values_merge.pt")

            self.gene_exp_graph_pre_exp_values_merge = gene_exp_graph_pre_exp_values_merge
        else:
            self.gene_exp_graph_pre_exp_values_merge = torch.load(f"./data/datastore/gene_exp_graph_pre_exp_values_merge.pt")

        for idx, sample in enumerate(self.dataset.iterrows()):
            if (not osp.exists(f"./data/datastore/data_{idx}{self.file_extension}.pt")):
                data = Data()
                data.x = torch.tensor(sample[1].values.astype(float)).reshape(-1,1).float()
                data.y = torch.tensor(self.labels.loc[sample[0], 'longevity'])
                data.age = torch.tensor(self.age.loc[sample[0], 'age'])

                if self.pre_filter is not None and not self.pre_filter(data):
                    continue

                if self.pre_transform is not None:
                    data = self.pre_transform(data)

                torch.save(data, osp.join(self.processed_dir, f'data_{idx}_{self.file_extension}.pt'))

    @property
    def processed_file_names(self):
        return [f"data_{i}{self.file_extension}.pt" for i in range(self.dataset.shape[0])]

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, f'data_{idx}_{self.file_extension}.pt'))
        return data