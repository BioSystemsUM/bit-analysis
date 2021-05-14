import os
from collections import defaultdict
from typing import List

import pandas as pd
from cobra import Model
from cobra.io import read_sbml_model
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler


def read_models(workdir: str) -> List[Model]:
    if not os.path.exists(workdir):
        raise OSError(f'{workdir} does not exist')

    models = os.listdir(workdir)

    cobra_models = []
    for model in models:

        if model.endswith('.xml') or model.endswith('.sbml'):
            model_path = os.path.join(workdir, model)
            cobra_model: Model = read_sbml_model(model_path)

            model_annotation = model.split('_')

            if len(model_annotation) > 1:
                _, organism, label = model_annotation
                cobra_model.id = f'{organism}_{label}'

            cobra_models.append(cobra_model)

    return cobra_models


def features_dataframe(models: List[Model],
                       factors: dict,
                       reactions: bool = True,
                       metabolites: bool = False) -> pd.DataFrame:
    if reactions and metabolites:
        raise ValueError('Select only one feature type, either reactions or metabolites')

    features_lookup = defaultdict(list)

    index = []
    for model in models:

        index.append(model.id)

        if reactions:
            for rxn in model.reactions:
                features_lookup[rxn.id].append(model.id)

        if metabolites:
            for met in model.metabolites:
                features_lookup[met.id].append(model.id)

    data = [[0] * len(features_lookup)] * len(index)
    df = pd.DataFrame(data=data,
                      index=index,
                      columns=features_lookup.keys())

    for key, models_ids in features_lookup.items():

        for model_id in models_ids:
            df.loc[model_id, key] = 1

    dfs = [df]
    for factor, values in factors.items():
        df = pd.DataFrame(data=values,
                          index=index,
                          columns=[factor])

        dfs.append(df)

    return pd.concat(dfs, axis=1)


def scaling(dataframe: pd.DataFrame, factors, standard=True, variance=True):
    x_mask = dataframe.columns[~dataframe.columns.isin(factors)]
    x = dataframe.loc[:, x_mask]
    y = dataframe.loc[:, factors]

    if variance:
        scalier = VarianceThreshold()
        scaled = scalier.fit(x)
        x = x.iloc[:, scaled.get_support(indices=True)]

    if standard:
        t_df = x.T

        scalier = StandardScaler()
        scaled = scalier.fit_transform(t_df)

        x = pd.DataFrame(scaled, columns=x.index, index=x.columns)
        x = x.T

    out_df = pd.concat([x, y], axis=1)

    return out_df


def pca_analysis(dataframe: pd.DataFrame, factors, components=2):
    x_mask = dataframe.columns[~dataframe.columns.isin(factors)]
    x = dataframe.loc[:, x_mask]
    y = dataframe.loc[:, factors]

    pca = PCA(n_components=components)

    pc = pca.fit_transform(x)

    columns = [f'PC {i + 1}' for i in range(components)]

    pc_df = pd.DataFrame(data=pc, index=dataframe.index, columns=columns)

    df = pd.concat([pc_df, y], axis=1)

    return df


def plot_pca(dataframe: pd.DataFrame, pc1, pc2, factors):
    for factor in factors:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(pc1, fontsize=15)
        ax.set_ylabel(pc2, fontsize=15)
        ax.set_title(f'{factor.title()} PCA', fontsize=20)
        targets = set(dataframe.loc[:, factor])

        for target, color in zip(targets, mcolors.TABLEAU_COLORS):
            idxs = dataframe.loc[:, factor] == target

            ax.scatter(dataframe.loc[idxs, pc1],
                       dataframe.loc[idxs, pc2],
                       s=50)

        ax.legend(targets)
        ax.grid()
        fig.show()
        fname = f'{factor.title()}_{pc1}_{pc2}_PCA.png'.replace(' ', '_')
        fig.savefig(fname=fname)


if __name__ == '__main__':
    directory = os.path.join(os.getcwd(), 'models')
    analysis_models = read_models(directory)
    models_factors = {'template': ['all', 'all', 'all'],
                      'organism': ['mtb', 'ec', 'ec']}
    models_factors_keys = tuple(models_factors.keys())
    rxns_df = features_dataframe(analysis_models, factors=models_factors, reactions=True)
    scaled_df = scaling(rxns_df, models_factors_keys)
    pca = pca_analysis(scaled_df, models_factors_keys)
    plot_pca(pca, 'PC 1', 'PC 2', models_factors_keys)
