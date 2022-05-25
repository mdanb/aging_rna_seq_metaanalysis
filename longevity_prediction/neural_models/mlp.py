from torch.nn import Sequential as Seq, Linear as Lin
import torch.nn.functional as F
from torch.nn import BatchNorm1d, ReLU, Dropout

class MLP(Seq):
    def __init__(self, input_dim, hidden_dim, output_dim, num_layers, dropout):
        layers = []
        layers.append(Lin(input_dim, hidden_dim))

        for i in range(num_layers - 2):
            layers.append(BatchNorm1d(hidden_dim))
            layers.append(ReLU())
            if dropout > 0:
                layers.append(Dropout(dropout))
            if (not (i == num_layers - 3)):
                layers.append(Lin(hidden_dim, hidden_dim))

        layers.append(Lin(hidden_dim, output_dim))

        self.layers = layers
        super(MLP, self).__init__(*self.layers)

    def loss(self, pred, label):
        return F.nll_loss(pred, label)

    def weight_reset(self, m):
        if isinstance(m, Lin):
            m.reset_parameters()

    def forward(self, x):
        for module in self:
            x = module(x)
        return x