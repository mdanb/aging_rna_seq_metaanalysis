import argparse
import torch
from utils import load_data, load_train_test_datasets, train_epoch, test_epoch, get_num_genes
from neural_models.mlp import MLP
from torch_geometric.data import DataLoader, Data

parser = argparse.ArgumentParser()

parser.add_argument('--expression_path', type=str,
                    default="../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv",
                    help='path to gene expression data '
                         '(default: ../common_datastore/getmm_combat_seq_no_outliers_and_singles_gene_expression.csv)')
parser.add_argument('--label_path', type=str, default="../common_datastore/labels.csv",
                    help='path to labels (default: ../common_datastore/labels.csv)')
parser.add_argument('--age_path', type=str, default="../common_datastore/age.csv",
                    help='path to age data (default: ../common_datastore/age.csv)')
parser.add_argument('--experiments_path', type=str, default="../common_datastore/sra_to_bioproject.csv",
                    help='path to sra to bioproject mapping (default: ../common_datastore/sra_to_bioproject.csv)')
# MLP parameters
parser.add_argument('--mlp_hidden_dim', type=int, default=512,
                    help='embedding dimensions (default: 512)')
parser.add_argument('--num_mlp_layers', type=int, default=3,
                    help='number of MLP layers total, excluding input layer (default: 3)')
# Training / Testing settings
parser.add_argument('--dropout', type=float, default=0,
                    help='dropout (default: 0)')
parser.add_argument('--weight_decay', type=float, default=0.001,
                    help='weight decay (default: 0.001)')
parser.add_argument('--batch_size', type=int, default=1500,
                    help='batch size for training (default: 1500)')
parser.add_argument('--learning_rate', type=float, default=0.0001,
                    help='learning rate (default: 0.0001)')
parser.add_argument('--epochs', type=int, default=100,
                    help='num training epochs (default: 100)')
parser.add_argument('--seed', type=int, default=42, help="Seed")
# Data Filtering
parser.add_argument('--aging_genes_only', action='store_true',
                    help="train the model using aging genes only (default: False)")

config = parser.parse_args()

epochs = config.epochs
batch_size = config.batch_size
learning_rate = config.learning_rate
weight_decay = config.weight_decay
seed = config.seed
mlp_hidden_dim = config.mlp_hidden_dim
data = config.expression_path.split('/')[-1]
aging_genes_only = config.aging_genes_only
num_mlp_layers = config.num_mlp_layers

device = torch.device("cuda")
torch.manual_seed(config.seed)

dataset = load_data(config)
train_test_dataset_list = load_train_test_datasets(dataset)
X_train, labels_train, _ = train_test_dataset_list[0]
X_test, labels_test, _ = train_test_dataset_list[1]
num_genes = get_num_genes(dataset)

MLP = MLP(num_genes + 1, config.mlp_hidden_dim, 3, config.num_mlp_layers,
          dropout=config.dropout).cuda()
opt = torch.optim.Adam(MLP.parameters(), lr=learning_rate, weight_decay=weight_decay)

xy_train = [Data(x=torch.Tensor(x.values), y=labels_train.loc[idx].values) for idx, x in X_train.iterrows()]
xy_test = [Data(x=torch.Tensor(x.values), y=labels_test.loc[idx].values) for idx, x in X_test.iterrows()]

train_loader = DataLoader(xy_train, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(xy_test, batch_size=batch_size, shuffle=True)

for epoch in range(1, epochs + 1):
    MLP, _, _ = train_epoch(train_loader, MLP, opt, device)

test_acc, _ = test_epoch(test_loader, MLP, device)

with open("results/neural/best_model_test_set_performance.txt", "w") as f:
    f.write(f"{test_acc}\n")