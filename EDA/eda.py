import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
import os
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.cm as mplcm
import math
from matplotlib.colors import ListedColormap
from matplotlib.cm import hsv
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

def pca(gene_exp_datasets, experiments, age, labels, subplots_dims, figsize, savefile, by_experiment=True):
    import matplotlib.colors as colors
    sp_dim1 = subplots_dims[0]
    sp_dim2 = subplots_dims[1]
    pca = PCA(n_components = 2)
    fig, ax = plt.subplots(sp_dim1, sp_dim2, figsize=(figsize[0], figsize[1]))

    for idx, (gene_expression_data, experiments, age, labels) in enumerate(zip(gene_exp_datasets, experiments, age, labels)):
        df_2d = pd.DataFrame(pca.fit_transform(gene_expression_data))
        df_2d.index = gene_expression_data.index

        if (by_experiment):
            # Taken from https://stackoverflow.com/questions/42697933/colormap-with-maximum-distinguishable-colours
            def generate_colormap(number_of_distinct_colors: int = 80):
                if number_of_distinct_colors == 0:
                    number_of_distinct_colors = 80

                number_of_shades = 7
                number_of_distinct_colors_with_multiply_of_shades = int(math.ceil(number_of_distinct_colors / number_of_shades) * number_of_shades)

                # Create an array with uniformly drawn floats taken from <0, 1) partition
                linearly_distributed_nums = np.arange(number_of_distinct_colors_with_multiply_of_shades) / number_of_distinct_colors_with_multiply_of_shades

                # We are going to reorganise monotonically growing numbers in such way that there will be single array with saw-like pattern
                #     but each saw tooth is slightly higher than the one before
                # First divide linearly_distributed_nums into number_of_shades sub-arrays containing linearly distributed numbers
                arr_by_shade_rows = linearly_distributed_nums.reshape(number_of_shades, number_of_distinct_colors_with_multiply_of_shades // number_of_shades)

                # Transpose the above matrix (columns become rows) - as a result each row contains saw tooth with values slightly higher than row above
                arr_by_shade_columns = arr_by_shade_rows.T

                # Keep number of saw teeth for later
                number_of_partitions = arr_by_shade_columns.shape[0]

                # Flatten the above matrix - join each row into single array
                nums_distributed_like_rising_saw = arr_by_shade_columns.reshape(-1)

                # HSV colour map is cyclic (https://matplotlib.org/tutorials/colors/colormaps.html#cyclic), we'll use this property
                initial_cm = hsv(nums_distributed_like_rising_saw)

                lower_partitions_half = number_of_partitions // 2
                upper_partitions_half = number_of_partitions - lower_partitions_half

                # Modify lower half in such way that colours towards beginning of partition are darker
                # First colours are affected more, colours closer to the middle are affected less
                lower_half = lower_partitions_half * number_of_shades
                for i in range(3):
                    initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8/lower_half)

                # Modify second half in such way that colours towards end of partition are less intense and brighter
                # Colours closer to the middle are affected less, colours closer to the end are affected more
                for i in range(3):
                    for j in range(upper_partitions_half):
                        modifier = np.ones(number_of_shades) - initial_cm[lower_half + j * number_of_shades: lower_half + (j + 1) * number_of_shades, i]
                        modifier = j * modifier / upper_partitions_half
                        initial_cm[lower_half + j * number_of_shades: lower_half + (j + 1) * number_of_shades, i] += modifier

                return ListedColormap(initial_cm)

            data_with_experiments = gene_expression_data.merge(experiments, left_index=True, right_index=True, how="left")
            num_colors = len(pd.unique(data_with_experiments.iloc[:, -1]))
            cm = generate_colormap(num_colors)
            cNorm  = colors.Normalize(vmin=0, vmax=num_colors)
            scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
            if len(gene_exp_datasets) > 1:
                if (sp_dim1 > 1):
                    ax[int(idx / sp_dim1)][idx % sp_dim2].set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colors)])
                else:
                    ax[idx].set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colors)])
            else:
                ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colors)])
            for _, group in data_with_experiments.groupby('Bioproject'):
                if len(gene_exp_datasets) > 1:
                    if (sp_dim1 > 1):
                        ax[int(idx / sp_dim1)][idx % sp_dim2].scatter(df_2d.loc[group.index][0], df_2d.loc[group.index][1], s=150)
                    else:
                        ax[idx].scatter(df_2d.loc[group.index][0], df_2d.loc[group.index][1], s=150)
                else:
                    ax.scatter(df_2d.loc[group.index][0], df_2d.loc[group.index][1], s=150)
        else:
            df_2d = df_2d.merge(age, left_index=True, right_index=True, how='left')
            df_2d = df_2d.merge(labels, left_index=True, right_index=True, how='left')
            df_2d.columns = ['dim1', 'dim2', 'age', 'longevity']
            cm = plt.get_cmap('plasma')
            data_with_age = gene_expression_data.merge(age, left_index=True, right_index=True, how="left")
            data_with_age = data_with_age.merge(labels, left_index=True, right_index=True, how="left")
            age_groups = df_2d.groupby('age')
            ages = list(np.sort(age['age'].unique())[::-1])
            colors = cm(np.linspace(0, 1, len(ages)))
            color_index = 0
            while (ages):
                for name, group in age_groups:
                    if(name == ages[0]):
                        if (len(gene_exp_datasets) > 1):
                            if (sp_dim1 > 1):
                                ax[int(idx / sp_dim1)][idx % sp_dim2].plot(group[group['longevity'] == 0].dim1, group[group['longevity'] == 0].dim2,
                                        label='long', marker='o', linestyle='',
                                       markersize=12, c=colors[color_index])
                                ax[int(idx / sp_dim1)][idx % sp_dim2].plot(group[group['longevity'] == 1].dim1, group[group['longevity'] == 1].dim2,
                                        label='normal', marker='*', linestyle='',
                                       markersize=12, c=colors[color_index])
                                ax[int(idx / sp_dim1)][idx % sp_dim2].plot(group[group['longevity'] == 2].dim1, group[group['longevity'] == 2].dim2,
                                        label='short', marker='^', linestyle='',
                                       markersize=12, c=colors[color_index])
                            else:
                                ax[idx].plot(group[group['longevity'] == 0].dim1, group[group['longevity'] == 0].dim2,
                                    label='long', marker='o', linestyle='',
                                   markersize=12, c=colors[color_index])
                                ax[idx].plot(group[group['longevity'] == 1].dim1, group[group['longevity'] == 1].dim2,
                                    label='normal', marker='*', linestyle='',
                                   markersize=12, c=colors[color_index])
                                ax[idx].plot(group[group['longevity'] == 2].dim1, group[group['longevity'] == 2].dim2,
                                    label='short', marker='^', linestyle='',
                                    markersize=12, c=colors[color_index])
                        else:
                            ax.plot(group[group['longevity'] == 0].dim1, group[group['longevity'] == 0].dim2,
                                    label='long', marker='o', linestyle='',
                                   markersize=12, c=colors[color_index])
                            ax.plot(group[group['longevity'] == 1].dim1, group[group['longevity'] == 1].dim2,
                                    label='normal', marker='*', linestyle='',
                                   markersize=12, c=colors[color_index])
                            ax.plot(group[group['longevity'] == 2].dim1, group[group['longevity'] == 2].dim2,
                                    label='short', marker='^', linestyle='',
                                   markersize=12, c=colors[color_index])

                        del ages[0]
                        color_index += 1
                        break

            ages = [str(age) + ' y/o' for age in list(np.sort(age['age'].unique())[::-1])]
            legend_elements = [Line2D([0], [0], marker='o', color='w', label='long',
                                      markerfacecolor='k', markersize=12),
                               Line2D([0], [0], marker='*', color='w', label='normal',
                                      markerfacecolor='k', markersize=12),
                               Line2D([0], [0], marker='^', color='w', label='short',
                                      markerfacecolor='k', markersize=12)]
            for i, color in enumerate(colors):
                age_patch = mpatches.Patch(color=colors[i], label=ages[i])
                legend_elements.append(age_patch)

            if (len(gene_exp_datasets) > 1):
                if (sp_dim1 > 1):
                    ax[int(idx / sp_dim1)][idx % sp_dim2].legend(handles=legend_elements)
                else:
                    ax[idx].legend(handles=legend_elements)
            else:
                ax.legend(handles=legend_elements)
    fig.savefig(f'figures/{savefile}.png', bbox_inches='tight')


def create_similarity_matrix(gene_expression_data, savefile, subplots_dims, figsize):
    fig, ax = plt.subplots(subplots_dims[0], subplots_dims[1], figsize=(figsize[0], figsize[1]))
    sp_dim1 = subplots_dims[0]
    sp_dim2 = subplots_dims[1]
    for idx, gene_expression in enumerate(gene_expression_data):
        distance_matrix = euclidean_distances(gene_expression)
        distance_matrix = pd.DataFrame(distance_matrix)
        np.fill_diagonal(distance_matrix.values, distance_matrix.max().max())
        distance_matrix = -distance_matrix
        if (len(gene_expression_data) > 1):
            if (sp_dim1 > 1):
                ax[int(idx / sp_dim1)][idx % sp_dim2].imshow(distance_matrix, cmap='hot')
            else:
                ax[idx].imshow(distance_matrix, cmap='hot')
        else:
            ax.imshow(distance_matrix, cmap='hot')

    fig.savefig(f'figures/{savefile}.png', bbox_inches='tight')
#     plt.imsave(f'figures/{savefile}.png', -distance_matrix, cmap='hot')
#         df_to_excel = (-distance_matrix).style.background_gradient(cmap='hot')
#     if not (os.path.exists(f'figures/{savefile}.xlsx')):
#         df_to_excel.to_excel(f'figures/{savefile}.xlsx', engine='openpyxl')

outliers = ['SRR10222650', 'SRR10222651', 'SRR10222652',
            'SRR10222653', 'SRR10222654', 'SRR10222655',
            'SRR10222656', 'SRR10222657', 'SRR10222658',
            'SRR10222659', 'SRR10222660', 'SRR10222661']

os.makedirs('figures', exist_ok=True)

experiments = pd.read_csv('../common_datastore/sra_to_bioproject.csv')
experiments = experiments.set_index("Sample")
age = pd.read_csv("../common_datastore/age.csv")
age.set_index('Sample', inplace=True)
labels = pd.read_csv("../common_datastore/labels.csv")
labels.set_index('Sample', inplace=True)

stats = labels[~labels.index.isin(outliers)].value_counts()
stats / sum(stats)

saucie_getmm_ge_no_outliers = pd.read_csv('../common_datastore/SAUCIE_getmm.csv')
saucie_getmm_ge_no_outliers.columns = ['Sample'] + list(saucie_getmm_ge_no_outliers.columns[1:])
saucie_getmm_ge_no_outliers = saucie_getmm_ge_no_outliers.set_index('Sample')

create_similarity_matrix([saucie_getmm_ge_no_outliers], "saucie_getmm_no_outliers_similarity_matrix", (1,1), (10,10))