import os
from collections import defaultdict
from typing import List, Tuple

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

from utils import (ModelAnalysis, parse_organism_id, read_models, parse_reaction, parse_metabolite,
                   parse_organism_annotation, parse_template_annotation, parse_method_annotation, COLORS,
                   get_explained_variance_idx)


def models_dataframe(models: List[ModelAnalysis],
                     reactions: bool = True,
                     filter_boundaries: bool = True) -> Tuple[pd.DataFrame, Tuple[str]]:
    features_lookup = defaultdict(list)

    categorical = defaultdict(list)
    index = []

    for model_analysis in models:

        index.append(model_analysis.model_id)
        categorical['organism_id'].append(model_analysis.organism_id)
        categorical['organism'].append(model_analysis.organism)
        categorical['template'].append(model_analysis.template)
        categorical['method'].append(model_analysis.method)

        if reactions:

            for rxn in model_analysis.model.reactions:

                rxn_id = parse_reaction(rxn, filter_boundaries)

                if rxn_id is not None:
                    features_lookup[rxn_id].append(model_analysis.model_id)

        else:

            for met in model_analysis.model.metabolites:

                met_id = parse_metabolite(met, filter_boundaries)

                if met_id is not None:
                    features_lookup[met_id].append(model_analysis.model_id)

    data = [[0] * len(features_lookup)] * len(index)
    df = pd.DataFrame(data=data,
                      index=index,
                      columns=features_lookup.keys())

    for rxn_or_met, models_lkp in features_lookup.items():

        for model_id in models_lkp:
            df.loc[model_id, rxn_or_met] = 1

    categorical_df = pd.DataFrame.from_dict(categorical)
    categorical_df.index = index

    df = pd.concat([df, categorical_df], axis=1)

    return df, tuple(categorical.keys())


def cog_dataframe(file_path: str) -> Tuple[pd.DataFrame, Tuple[str]]:
    df = pd.read_csv(file_path,
                     sep='\t',
                     index_col='organism')
    df.loc[:, 'organism_id'] = [parse_organism_id(organism) for organism in df.index]
    return df, ('domain', 'phylum', 'organism_id')


def models_genes_cog_dataframe(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path,
                     sep='\t',
                     index_col='model_id')

    organisms = []
    organisms_id = []
    templates = []
    methods = []
    for model_id in df.index:
        model_annotation = model_id.split('_')

        organism = parse_organism_annotation(model_annotation)
        organism_id = parse_organism_id(organism)
        template = parse_template_annotation(model_annotation)
        method = parse_method_annotation(model_annotation)

        organisms.append(organism)
        organisms_id.append(organism_id)
        templates.append(template)
        methods.append(method)

    df.loc[:, 'organism_id'] = organisms_id
    df.loc[:, 'organism'] = organisms
    df.loc[:, 'template'] = templates
    df.loc[:, 'method'] = methods

    return df, ('organism_id', 'organism', 'template', 'method')


def scaling(dataframe: pd.DataFrame,
            categorical: Tuple[str],
            standard: bool = True,
            variance: bool = True) -> pd.DataFrame:
    mask = dataframe.columns.isin(categorical)
    cols = dataframe.columns[~mask]
    x = dataframe.loc[:, cols]
    y = dataframe.loc[:, categorical]

    if variance:
        scalier = VarianceThreshold()
        scaled = scalier.fit(x)
        x = x.iloc[:, scaled.get_support(indices=True)]

    if standard:
        scalier = StandardScaler()
        scaled = scalier.fit_transform(x.T)

        x = pd.DataFrame(scaled, columns=x.index, index=x.columns)
        x = x.T

    return pd.concat([x, y], axis=1)


def method_filter(dataframe: pd.DataFrame, label: str):
    mask = dataframe.loc[:, 'method'] != label

    return dataframe.loc[mask, :]


def pca_analysis(dataframe: pd.DataFrame,
                 categorical: Tuple[str],
                 components: int = 2) -> Tuple[pd.DataFrame, PCA]:
    mask = dataframe.columns.isin(categorical)
    cols = dataframe.columns[~mask]
    x = dataframe.loc[:, cols]
    y = dataframe.loc[:, categorical]

    pca = PCA(n_components=components)

    pc = pca.fit_transform(x)

    columns = [f'PC {i + 1}' for i in range(components)]

    df = pd.DataFrame(data=pc, index=dataframe.index, columns=columns)

    df = pd.concat([df, y], axis=1)
    return df, pca


def plot_pca(workdir: str,
             dataframe: pd.DataFrame,
             pca: PCA,
             c1: str,
             c2: str,
             title: str,
             factor: str):
    explained_variance_1 = get_explained_variance_idx(c1)
    explained_variance_2 = get_explained_variance_idx(c2)

    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title, fontsize=20)

    x_label = f'{c1} ({round(pca.explained_variance_ratio_[explained_variance_1] * 100, 2)} %)'
    y_label = f'{c2} ({round(pca.explained_variance_ratio_[explained_variance_2] * 100, 2)} %)'
    ax.set_xlabel(x_label, fontsize=15)
    ax.set_ylabel(y_label, fontsize=15)

    labels = set(dataframe.loc[:, factor])

    for label, color in zip(labels, COLORS):
        mask = dataframe.loc[:, factor] == label

        pc1 = dataframe.loc[mask, c1]
        pc2 = dataframe.loc[mask, c2]
        ax.scatter(pc1,
                   pc2,
                   c=color,
                   s=60)

        organisms_id = dataframe.loc[mask, 'organism_id']

        for pc1_pt, pc2_pt, annotation in zip(pc1, pc2, organisms_id):
            ax.annotate(annotation, (pc1_pt + 1, pc2_pt - 1))

    legend = ax.legend(labels, loc=(1.04, 0))
    ax.grid()
    file_name = f'{title}_{factor}_{c1}_{c2}.png'
    file_path = os.path.join(workdir, file_name)
    fig.savefig(fname=file_path, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)


# ----------------------------------
# ANALYSIS RUN
# ----------------------------------
def organisms_analysis(cog_file: str,
                       analysis_dir: str):
    df, categorical = cog_dataframe(cog_file)
    df = scaling(dataframe=df, categorical=categorical)
    pca_df, pca = pca_analysis(dataframe=df, categorical=categorical, components=3)

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='Metabolic COG Analysis',
             factor='domain')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='Metabolic COG Analysis',
             factor='domain')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='Metabolic COG Analysis',
             factor='phylum')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='Metabolic COG Analysis',
             factor='phylum')


def reactions_analysis(models_dir: str,
                       analysis_dir: str,
                       filter_boundaries: bool,
                       method: str = 'permissive',
                       read: str = '',
                       write: str = ''):
    if read:

        df = pd.read_csv(read,
                         sep='\t',
                         index_col='model_id')

        categorical = ('organism_id', 'organism', 'template', 'method')

    else:
        models_analysis = read_models(models_dir)
        df, categorical = models_dataframe(models_analysis, reactions=True, filter_boundaries=filter_boundaries)

        if write:
            df.to_csv(write, sep='\t', index=True, index_label='model_id')

    df = method_filter(df, method)

    df = scaling(dataframe=df, categorical=categorical)
    pca_df, pca = pca_analysis(dataframe=df, categorical=categorical, components=3)

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='Reactions Analysis',
             factor='template')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='Reactions Analysis',
             factor='template')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='Reactions Analysis',
             factor='method')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='Reactions Analysis',
             factor='method')


def metabolites_analysis(models_dir: str,
                         analysis_dir: str,
                         filter_boundaries: bool,
                         method: str = 'permissive',
                         read: str = '',
                         write: str = ''):
    if read:
        df = pd.read_csv(read,
                         sep='\t',
                         index_col='model_id')

        categorical = ('organism_id', 'organism', 'template', 'method')

    else:
        models_analysis = read_models(models_dir)
        df, categorical = models_dataframe(models_analysis, reactions=False, filter_boundaries=filter_boundaries)

        if write:
            df.to_csv(write, sep='\t', index=True, index_label='model_id')

    df = method_filter(df, method)

    df = scaling(dataframe=df, categorical=categorical)
    pca_df, pca = pca_analysis(dataframe=df, categorical=categorical, components=3)

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='Metabolites Analysis',
             factor='template')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='Metabolites Analysis',
             factor='template')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='Metabolites Analysis',
             factor='method')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='Metabolites Analysis',
             factor='method')


def genes_analysis(cog_file: str, analysis_dir: str, method: str = 'permissive'):

    df, categorical = models_genes_cog_dataframe(cog_file)

    df = method_filter(df, method)

    df = scaling(df, categorical=categorical)
    pca_df, pca = pca_analysis(df, categorical=categorical, components=3)

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='COG Genes Analysis',
             factor='template')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='COG Genes Analysis',
             factor='template')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 2',
             title='COG Genes Analysis',
             factor='method')

    plot_pca(workdir=analysis_dir,
             dataframe=pca_df,
             pca=pca,
             c1='PC 1',
             c2='PC 3',
             title='COG Genes Analysis',
             factor='method')


if __name__ == '__main__':
    base_dir = os.getcwd()
    models_base_dir = os.path.join(base_dir, 'models')
    comparative_dir = os.path.join(base_dir, 'comparative_analysis')
    models_analysis_dir = os.path.join(base_dir, 'model_analysis')

    organisms_cog_file = os.path.join(comparative_dir, 'genomes_cog_analysis.tsv')
    models_genes_cog_file = os.path.join(comparative_dir, 'models_cog_analysis.tsv')

    reactions_dir = os.path.join(models_analysis_dir, 'reactions.tsv')
    metabolites_dir = os.path.join(models_analysis_dir, 'metabolites.tsv')

    organisms_analysis(cog_file=organisms_cog_file,
                       analysis_dir=comparative_dir)

    reactions_analysis(models_dir=models_base_dir,
                       analysis_dir=models_analysis_dir,
                       filter_boundaries=True,
                       write=reactions_dir)

    metabolites_analysis(models_dir=models_base_dir,
                         analysis_dir=models_analysis_dir,
                         filter_boundaries=True,
                         write=metabolites_dir)

    # genes_analysis(cog_file=models_genes_cog_file,
    #                analysis_dir=comparative_dir)
