import os
from collections import defaultdict
from typing import List, Union, Tuple

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

from utils import ModelAnalysis, parse_organism_id, read_models, parse_reaction, parse_metabolite, \
    parse_organism_annotation, parse_template_annotation, parse_method_annotation, COLORS


def models_dataframe(models: List[ModelAnalysis],
                     reactions: bool = True,
                     filter_boundaries: bool = True,
                     write: str = '') -> pd.DataFrame:
    features_lookup = defaultdict(list)

    models_ids = []
    organisms = []
    organisms_id = []
    templates = []
    methods = []
    for model_analysis in models:

        models_ids.append(model_analysis.model_id)
        organisms.append(model_analysis.organism)
        organisms_id.append(model_analysis.organism_id)
        templates.append(model_analysis.template)
        methods.append(model_analysis.method)

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

    data = [[0] * len(features_lookup)] * len(models_ids)
    df = pd.DataFrame(data=data,
                      index=models_ids,
                      columns=features_lookup.keys())

    for rxn_or_met, models_lkp in features_lookup.items():

        for model_id in models_lkp:
            df.loc[model_id, rxn_or_met] = 1

    dfs = [df]
    factors = ('model_id', 'organism', 'organism_id', 'template', 'method')
    labels = (models_ids, organisms, organisms_id, templates, methods)

    for factor, label in zip(factors, labels):
        df = pd.DataFrame(data=label,
                          index=models_ids,
                          columns=[factor])

        dfs.append(df)

    df = pd.concat(dfs, axis=1)

    if write:
        if reactions:
            df.to_excel(write,
                        sheet_name='reactions',
                        index_label='id')

        else:
            df.to_excel(write,
                        sheet_name='metabolites',
                        index_label='id')

    return df


def cog_dataframe(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path,
                     sep='\t',
                     index_col='organism')
    df.loc[:, 'organism_id'] = [parse_organism_id(organism) for organism in df.index]
    return df


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

    df.loc[:, 'organism'] = organisms
    df.loc[:, 'organism_id'] = organisms_id
    df.loc[:, 'template'] = templates
    df.loc[:, 'method'] = methods

    return df


def scaling(dataframe: pd.DataFrame,
            factors: Union[List[str], Tuple[str]],
            model_filter: float = 0.0,
            standard: bool = True,
            variance: bool = True) -> pd.DataFrame:

    # filter rows
    dataframe = dataframe[dataframe.sum(axis=1) > round(dataframe.shape[1] * model_filter, 0)]

    x_mask = dataframe.columns[~dataframe.columns.isin(factors)]
    x = dataframe.loc[:, x_mask]
    y = dataframe.loc[:, factors]

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


def pca_analysis(dataframe: pd.DataFrame,
                 factors: Union[List[str], Tuple[str]],
                 components: int = 2) -> Tuple[pd.DataFrame, PCA]:

    x_mask = dataframe.columns[~dataframe.columns.isin(factors)]
    x = dataframe.loc[:, x_mask]
    y = dataframe.loc[:, factors]

    pca = PCA(n_components=components)

    pc = pca.fit_transform(x)

    columns = [f'PC {i + 1}' for i in range(components)]

    df = pd.DataFrame(data=pc, index=dataframe.index, columns=columns)

    df = pd.concat([df, y], axis=1)
    return df, pca


def get_explained_variance_idx(component):
    if '1' in component:
        return 0

    elif '2' in component:
        return 1

    elif '3' in component:
        return 2

    return 3


def plot_pca(workdir: str,
             pca: Tuple[pd.DataFrame, PCA],
             c1: str,
             c2: str,
             factors: Union[List[str], Tuple[str]],
             content: str):

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    df, pca = pca

    explained_variance_1 = get_explained_variance_idx(c1)
    explained_variance_2 = get_explained_variance_idx(c2)

    for factor in factors:

        fig = plt.figure(figsize=(8, 8))

        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(content, fontsize=20)

        x_label = f'{c1} ({round(pca.explained_variance_ratio_[explained_variance_1] * 100, 2)} %)'
        y_label = f'{c2} ({round(pca.explained_variance_ratio_[explained_variance_2] * 100, 2)} %)'
        ax.set_xlabel(x_label, fontsize=15)
        ax.set_ylabel(y_label, fontsize=15)

        labels = set(df.loc[:, factor])

        for label, color in zip(labels, COLORS):
            mask = df.loc[:, factor] == label

            pc1 = df.loc[mask, c1]
            pc2 = df.loc[mask, c2]
            ax.scatter(pc1,
                       pc2,
                       c=color,
                       s=60)

            organisms_id = df.loc[mask, 'organism_id']

            for pc1_pt, pc2_pt, annotation in zip(pc1, pc2, organisms_id):
                ax.annotate(annotation, (pc1_pt + 1, pc2_pt - 1))

        legend = ax.legend(labels, loc=(1.04, 0))
        ax.grid()
        file_name = f'{content}_{factor.title()}_{c1}_{c2}.png'
        file_path = os.path.join(workdir, file_name)
        fig.savefig(fname=file_path, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)


# ----------------------------------
# ANALYSIS RUN
# ----------------------------------
def organisms_func_analysis(cog_analysis_file: str,
                            analysis_dir: str):

    factors = ('domain', 'phylum', 'organism_id')

    df = cog_dataframe(cog_analysis_file)
    df = scaling(df, factors)
    pca = pca_analysis(dataframe=df, factors=factors, components=3)

    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 2',
             factors=('phylum', 'domain'), content='Metabolic COG Analysis')
    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 3',
             factors=('phylum', 'domain'), content='Metabolic COG Analysis')


def reactions_analysis(models_dir: str,
                       analysis_dir: str,
                       filter_exchanges: bool,
                       read: str = '',
                       write: str = ''):

    factors = ('model_id', 'organism', 'organism_id', 'template', 'method')

    if read:
        df = pd.read_excel(read, sheet_name='reactions', index_col='id')

    else:
        analysis_models = read_models(models_dir)
        df = models_dataframe(analysis_models, reactions=True, filter_boundaries=filter_exchanges, write=write)

    df = scaling(dataframe=df, factors=factors)
    pca = pca_analysis(dataframe=df, factors=factors, components=3)

    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 2',
             factors=('template', 'method'), content='Reactions Analysis')

    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 3',
             factors=('template', 'method'), content='Reactions Analysis')


def metabolites_analysis(models_dir: str,
                         analysis_dir: str,
                         filter_exchanges: bool,
                         read: str = '',
                         write: str = ''):
    factors = ('model_id', 'organism', 'organism_id', 'template', 'method')

    if read:
        df = pd.read_excel(read, sheet_name='metabolites', index_col='id')

    else:
        analysis_models = read_models(models_dir)
        df = models_dataframe(analysis_models, reactions=False, filter_boundaries=filter_exchanges, write=write)

    df = scaling(df, factors=factors)
    pca = pca_analysis(df, factors=factors, components=3)

    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 2',
             factors=('template', 'method'), content='Metabolites Analysis')

    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 3',
             factors=('template', 'method'), content='Metabolites Analysis')


def models_genes_cog_pca(cog_analysis_file: str, analysis_dir: str):
    factors = ('organism', 'organism_id', 'template', 'method')

    df = models_genes_cog_dataframe(cog_analysis_file)
    df = scaling(df, factors, 0)
    pca = pca_analysis(df, factors=factors, components=3)

    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 2',
             factors=('template', 'method'), content='Models Genes COG Analysis')
    plot_pca(workdir=analysis_dir, pca=pca, c1='PC 1', c2='PC 3',
             factors=('template', 'method'), content='Models Genes COG Analysis')


if __name__ == '__main__':

    # organisms_func_analysis(cog_analysis_file=os.path.join(os.getcwd(), 'comparative_func_analysis',
    #                                                        'genomes_cog_analysis.tsv'),
    #                         analysis_dir=os.path.join(os.getcwd(), 'comparative_func_analysis'))
    #
    # reactions_analysis(models_dir=os.path.join(os.getcwd(), 'models'),
    #                    analysis_dir=os.path.join(os.getcwd(), 'model_content_analysis'),
    #                    filter_exchanges=True,
    #                    read=os.path.join(os.getcwd(), 'model_content_analysis', 'reactions_analysis.xlsx'))
    #
    # metabolites_analysis(models_dir=os.path.join(os.getcwd(), 'models'),
    #                      analysis_dir=os.path.join(os.getcwd(), 'model_content_analysis'),
    #                      filter_exchanges=True,
    #                      read=os.path.join(os.getcwd(), 'model_content_analysis', 'metabolites_analysis.xlsx'), )

    models_genes_cog_pca(cog_analysis_file=os.path.join(os.getcwd(), 'comparative_func_analysis',
                                                        'models_cog_analysis.tsv'),
                         analysis_dir=os.path.join(os.getcwd(), 'model_content_analysis'))


