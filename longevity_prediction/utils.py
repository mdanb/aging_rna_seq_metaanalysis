import numpy as np
import torch
import os
from torch_geometric.data import DataLoader, Data
from tqdm import tqdm
from split import StratifiedGroupKFold
from pathlib import Path
from sklearn.model_selection import StratifiedShuffleSplit
import pandas as pd
from data.gene_expression_dataset import DatasetFromCSV
from data.gene_interaction_graph import StringDBGraph
from data.gene_expression_graph import GeneExpressionGraph

def cross_validate(model, X, labels, device, config, experiments, experiment_name,
                   num_folds, learning_rate, batch_size, weight_decay, epochs,
                   mixsplit=False):
    Path(f"gnn_checkpoints").mkdir(exist_ok=True)
    if (mixsplit):
        stratified_k_fold = StratifiedShuffleSplit(num_folds, test_size=0.1, random_state=config.seed)
    else:
        stratified_k_fold = StratifiedGroupKFold(num_folds)

    device = device
    metrics = {'all_val_accs': [], 'all_losses': []}
    fold_predictions = []
    for fold, (train_index, val_index) in enumerate(stratified_k_fold.split(np.zeros(len(labels)),
                                                                            labels.values.reshape(-1), experiments)):
        fold_data = pd.DataFrame(columns=['Sample', 'Experiment', 'Label', 'Prediction'])
        fold_data['Sample'] = experiments.index[val_index]
        fold_data['Experiment'] = experiments.values[val_index]
        fold_data['Label'] = labels['longevity']

        Path(f"./results/neural/{experiment_name}").mkdir(parents=True, exist_ok=True)
        f = open(f"./results/neural/{experiment_name}/fold_{fold}_training_progress.txt", "w")
        model.apply(model.weight_reset)
        if (config.train_GNN):
            xy_train_fold = X[list(train_index)]
            xy_val_fold = X[list(val_index)]
        else:
            xy_train_fold = [Data(x=torch.Tensor(x.values),
                                  y=labels.loc[idx].values) for idx, x in X.iloc[list(train_index)].iterrows()]
            xy_val_fold = [Data(x=torch.Tensor(x.values), y=labels.loc[idx].values) for idx, x in
                           X.iloc[list(val_index)].iterrows()]
        opt = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

        train_loader = DataLoader(xy_train_fold, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(xy_val_fold, batch_size=batch_size, shuffle=True)

        fold_val_accs = []
        fold_losses = []
        all_epoch_predictions = []
        loaded_upper_checkpoint = False
        for epoch in range(1, epochs + 1):
            upper_checkpoint = int(
                (epoch + config.eval_model_every) / config.eval_model_every) * config.eval_model_every \
                if epoch % config.eval_model_every != 0 and epoch != 1 else epoch
            upper_checkpoint_path = f"gnn_checkpoints/fold_{fold}_epoch_{upper_checkpoint}_{experiment_name}.pt.tar"
            path = f"gnn_checkpoints/fold_{fold}_epoch_{epoch}_{experiment_name}.pt.tar"

            if (os.path.exists(upper_checkpoint_path)):
                if (not loaded_upper_checkpoint):
                    checkpoint = torch.load(upper_checkpoint_path)
                    model.load_state_dict(checkpoint['model_state_dict'])
                    opt.load_state_dict(checkpoint['optimizer_state_dict'])
                    total_loss = checkpoint['loss']
                    loaded_upper_checkpoint = True
            else:
                model, opt, total_loss = train_epoch(train_loader, model, opt, device)

            if epoch % config.eval_model_every == 0 or epoch == 1:
                if (config.train_GNN):
                    if not loaded_upper_checkpoint:
                        torch.save({
                            'epoch': epoch,
                            'model_state_dict': model.state_dict(),
                            'optimizer_state_dict': opt.state_dict(),
                            'loss': total_loss,
                        }, path)

                    loaded_upper_checkpoint = False
                val_acc, epoch_predictions = test_epoch(val_loader, model, device)
                fold_val_accs.append(val_acc)
                fold_losses.append(total_loss)
                all_epoch_predictions.append(epoch_predictions)
                print(f"Fold {fold}. Epoch {epoch}. Loss: {total_loss}. Validation accuracy: {val_acc}")
                f.write(f"Fold {fold}. Epoch {epoch}. Loss: {total_loss}. Validation accuracy: {val_acc} \n")

        metrics['all_val_accs'].append(fold_val_accs)
        metrics['all_losses'].append(fold_losses)
        fold_predictions.append(all_epoch_predictions)

    avg_val_accs_across_folds = [float(sum(col)) / len(col) for col in zip(*metrics['all_val_accs'])]
    avg_loss_across_folds = [float(sum(col)) / len(col) for col in zip(*metrics['all_losses'])]

    max_val_acc = max(avg_val_accs_across_folds)
    max_val_index = avg_val_accs_across_folds.index(max_val_acc)
    loss = avg_loss_across_folds[max_val_index]
    best_epoch = max_val_index * 10 if max_val_index != 1 else 1

    individual_fold_performances = [row[int(best_epoch / config.eval_model_every)] for row in metrics['all_val_accs']]
    best_predictions_per_fold = [preds[int(best_epoch / config.eval_model_every)] for preds in fold_predictions]
    best_preds = [pred for preds in best_predictions_per_fold for pred in preds]
    best_preds = pd.DataFrame(best_preds)
    predicted_proportions = pd.DataFrame(best_preds.value_counts() / sum(best_preds.value_counts()),
                                         columns=['Predicted Proportions'])
    predicted_proportions.rename(index={0: "long-lived", 1: "normal-lived", 2: "short-lived"}, inplace=True)
    predicted_proportions.reset_index(inplace=True)
    predicted_proportions.columns = ['Longevity', 'Predicted Proportions']
    f = open(f"./results/neural/{experiment_name}/best_model_stats.txt", "a")
    print(f"Best Epoch {best_epoch}. Loss: {loss}. Best Average Validation Accuracy: {max_val_acc}")
    f.write(f"Best Epoch {best_epoch}. Loss: {loss}. Best Average Validation Accuracy: {max_val_acc} \n")
    f.write("\n")
    if (num_folds == 5):
        df_fold_performance = pd.DataFrame(np.zeros([1, 5]), columns=['Fold 0', 'Fold 1', 'Fold 2', 'Fold 3', 'Fold 4'])
    elif (num_folds == 10):
        df_fold_performance = pd.DataFrame(np.zeros([1, 10]), columns=['Fold 0', 'Fold 1', 'Fold 2', 'Fold 3', 'Fold 4',
                                                                       'Fold 5', 'Fold 6', 'Fold 7', 'Fold 8',
                                                                       'Fold 9'])
    df_fold_performance.iloc[0, :] = individual_fold_performances
    f.write(f"Performance across folds:\n")
    f.write(df_fold_performance.to_string(index=False))
    f.write('\n\n')
    f.write("Overall predicted proportions:\n")
    f.write(predicted_proportions.to_string(index=False))
    f.write('\n')

def train_epoch(loader, model, opt, device):
    total_loss = 0
    model.train()
    for batch in tqdm(loader):
        opt.zero_grad()
        label = batch.y
        if (model.__class__.__name__ == "MLP"):
            x = batch.x.reshape(len(batch.y), -1).to(device)
            pred = model(x)
            label = torch.LongTensor(np.concatenate(label))
        else:
            batch.to(device)
            pred = model(batch)

        loss = model.loss(pred, label.to(device))
        loss.backward()
        opt.step()
        total_loss += loss.item() * batch.num_graphs
    total_loss /= len(loader.dataset)
    return model, opt, total_loss

def test_epoch(loader, model, device):
    model.eval()
    correct = 0
    predictions = []
    for data in tqdm(loader):
        with torch.no_grad():
            label = data.y
            if (model.__class__.__name__ == "MLP"):
                x = data.x.reshape(len(data.y), -1).to(device)
                pred = model(x)
                label = torch.LongTensor(np.concatenate(label))
            else:
                data = data.to(device)
                pred = model(data)
            label = label.to(device)
            pred = pred.argmax(dim=1)
            predictions.extend(pred.tolist())

        correct += pred.eq(label).sum().item()

    total = len(loader.dataset)

    return correct / total, predictions

def load_data(config):
    dataset = DatasetFromCSV(f'{config.expression_path}',
                                 f'{config.label_path}',
                                 f'{config.age_path}',
                                 f'{config.experiments_path}',
                                 config.aging_genes_only)

    try:
        if (config.train_GNN):
            gene_interaction_graph = StringDBGraph()
            gene_expression_graphs = GeneExpressionGraph(gene_interaction_graph, dataset.df, config, dataset.labels,
                                                         dataset.age)
            pre_expression_merge_graph = gene_expression_graphs.gene_exp_graph_pre_exp_values_merge
            return dataset, gene_expression_graphs, pre_expression_merge_graph
        else:
            return dataset
    except:
        return dataset

def load_train_test_datasets(dataset, gene_expression_graphs=None):
    cv = StratifiedGroupKFold(10)
    labels = dataset.labels.to_numpy().reshape(-1)

    for idx, (train_idxs, test_idxs) in enumerate(cv.split(dataset.df.to_numpy(), labels,
                                                           dataset.experiments)):
        if (idx == 1):
            break
        if (gene_expression_graphs):
            X_train = gene_expression_graphs[train_idxs]
            X_test = gene_expression_graphs[test_idxs]
        else:
            age = dataset.age
            X = dataset.df.merge(age, left_index=True, right_index=True)
            X_train = X.iloc[train_idxs]
            X_test = X.iloc[test_idxs]

        labels_train = dataset.labels.iloc[train_idxs]
        experiments_train = dataset.experiments.iloc[train_idxs]
        labels_test = dataset.labels.iloc[test_idxs]
        experiments_test = dataset.experiments.iloc[test_idxs]

    return [(X_train, labels_train, experiments_train),
            (X_test, labels_test, experiments_test)]

def get_num_genes(dataset):
    return dataset.df.shape[1]