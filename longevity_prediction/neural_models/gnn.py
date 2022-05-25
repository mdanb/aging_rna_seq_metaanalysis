import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from torch_sparse import SparseTensor
from torch_geometric.nn import global_sort_pool
from torch.nn import BatchNorm1d
from neural_models.mlp import MLP

class GNN(torch.nn.Module):
    def __init__(self, input_dim, backbone_channels, output_dim, num_backbone_layers, concat_input_graph,
                 num_nodes, edge_attr, edge_index, k, mlp_hidden_dim, mlp_num_layers, dropout):
        super(GNN, self).__init__()
        edge_attr = torch.Tensor(edge_attr).float()
        self.input_dim = input_dim
        self.num_nodes = num_nodes
        self.edge_index_sp = SparseTensor(row=edge_index[0], col=edge_index[1],
                                          value=edge_attr, sparse_sizes=(num_nodes, num_nodes)).cuda()
        self.k = k
        self.backbone_channels = backbone_channels
        self.num_backbone_layers = num_backbone_layers
        self.dropout = dropout
        self.convs = torch.nn.ModuleList()
        self.convs.append(self.build_conv_model(input_dim, backbone_channels))
        self.bns = torch.nn.ModuleList()
        self.concat_input_graph = concat_input_graph
        for l in range(num_backbone_layers - 1):
            self.convs.append(self.build_conv_model(backbone_channels, backbone_channels))
            self.bns.append(BatchNorm1d(backbone_channels))

        if self.concat_input_graph:
            self.post_mp = MLP(k * num_backbone_layers * backbone_channels + k * input_dim + 1, mlp_hidden_dim, output_dim,
                               mlp_num_layers, dropout=dropout)
        else:
            self.post_mp = MLP(k * num_backbone_layers * backbone_channels + 1, mlp_hidden_dim, output_dim, mlp_num_layers,
                               dropout=dropout)

    def build_conv_model(self, input_dim, hidden_dim):
        return GCNConv(input_dim, hidden_dim)

    def forward(self, data):
        x, batch, edge_index, _ = data.x, data.batch, self.edge_index_sp, data.age
        batch_size = data.y.shape[0]
        x = x.reshape(batch_size, self.num_nodes, self.input_dim)

        feats_at_different_layers = []
        feats_at_different_layers.append(x)

        for i in range(self.num_backbone_layers):
            x = self.convs[i](x, edge_index)
            if not i == self.num_backbone_layers - 1:
                x = self.bns[i](x.reshape(-1, self.backbone_channels))
            x = x.reshape(batch_size, self.num_nodes, self.backbone_channels)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
            feats_at_different_layers.append(x)

        if self.concat_input_graph:
            x = torch.cat(tuple(feats_at_different_layers), dim=-1)
        else:
            x = torch.cat(tuple(feats_at_different_layers[1:]), dim=-1)

        if self.k != self.num_nodes:
            x = global_sort_pool(x.reshape(-1, self.num_backbone_layers * self.backbone_channels + self.input_dim),
                                 batch, self.k)
        if self.concat_input_graph:
            x = x.reshape(batch_size, self.k * (self.num_backbone_layers * self.backbone_channels + self.input_dim), 1)
        else:
            x = x.reshape(batch_size, self.k * (self.num_backbone_layers * self.backbone_channels), 1)
        data.x = x

        x = data.x.reshape(len(data.y), -1)
        age = data.age
        x = torch.cat((x, age.reshape(-1, 1)), dim=1)
        x = self.post_mp(x)

        return x

    def model_reset(self):
        for layer in self.children():
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()

    def weight_reset(self, m):
        if isinstance(m, torch.nn.Linear):
            m.reset_parameters()

    def loss(self, pred, label):
        return F.nll_loss(pred, label)