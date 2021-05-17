import os
from collections import defaultdict, namedtuple
from typing import List, Union, Tuple

import pandas as pd
from cobra import Model
from cobra.io import read_sbml_model
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

ModelAnalysis = namedtuple('ModelAnalysis',
                           field_names=('model',
                                        'organism',
                                        'organism_id',
                                        'template',
                                        'method'),
                           defaults=(Model(),
                                     'Mycobacterium tuberculosis',
                                     'Mtub'
                                     'all',
                                     'permissive'))


def parse_organism_annotation(model_annotation):
    if len(model_annotation) > 1:

        organism, *_ = model_annotation

        if 'tuber' in organism:
            return 'Mycobacterium tuberculosis'

        elif 'thermo' in organism:
            return 'Streptococcus thermophilus'

        elif 'fasti' in organism:
            return 'Xylella fastidiosa'

    return '_'.join(model_annotation)


def parse_organism_id(organism_name):
    genus_species = organism_name.split()
    genus = genus_species[0]
    species = genus_species[1]

    return f'{genus[0]}{species[0:3]}'


def parse_template_annotation(model_annotation):
    if len(model_annotation) > 1:

        _, template, *_ = model_annotation

        if 'all' in template:
            return 'all'

        elif 'random' in template:
            return 'random'

        elif 'select' in template:
            return 'select'

    return '_'.join(model_annotation)


def parse_method_annotation(model_annotation):
    if len(model_annotation) > 1:

        *_, method = model_annotation

        if 'permissive' in method:
            return 'permissive'

        elif 'restrictive' in method:
            return 'restrictive'

    return '_'.join(model_annotation)


def read_models(workdir: str) -> List[ModelAnalysis]:
    if not os.path.exists(workdir):
        raise OSError(f'{workdir} does not exist')

    models = os.listdir(workdir)

    models_analysis = []

    for model_name in models:

        if model_name.endswith('.xml') or model_name.endswith('.sbml'):
            model_path = os.path.join(workdir, model_name)

            model: Model = read_sbml_model(model_path)

            name = model_name.replace('model_', '')

            model_annotation = name.split('_')

            organism = parse_organism_annotation(model_annotation)
            organism_id = parse_organism_id(organism)
            template = parse_template_annotation(model_annotation)
            method = parse_method_annotation(model_annotation)

            model.id = f'{organism_id}_{template}_{method}'

            model_analysis = ModelAnalysis(model=model,
                                           organism=organism,
                                           organism_id=organism_id,
                                           template=template,
                                           method=method)

            models_analysis.append(model_analysis)

    return models_analysis


def models_dataframe(models: List[ModelAnalysis],
                     reactions: bool = True,
                     filter_exchanges: bool = True) -> pd.DataFrame:

    features_lookup = defaultdict(list)

    index = []
    organisms = []
    organisms_id = []
    templates = []
    methods = []
    for model_analysis in models:

        index.append(model_analysis.model.id)
        organisms.append(model_analysis.organism)
        organisms_id.append(model_analysis.organism_id)
        templates.append(model_analysis.template)
        methods.append(model_analysis.method)

        if reactions:

            for rxn in model_analysis.model.reactions:

                is_exchange = False

                if rxn.boundary:
                    is_exchange = True

                else:

                    for met in rxn.metabolites:
                        if met.id.endswith('_b'):
                            is_exchange = True
                            break

                if is_exchange:
                    continue

                else:
                    features_lookup[rxn.id].append(model_analysis.model.id)

        else:

            for met in model_analysis.model.metabolites:
                if met.id.endswith('_b'):
                    continue

                else:
                    features_lookup[met.id].append(model_analysis.model.id)

    data = [[0] * len(features_lookup)] * len(index)
    df = pd.DataFrame(data=data,
                      index=index,
                      columns=features_lookup.keys())

    for rxn_or_met, models_ids in features_lookup.items():

        for model_id in models_ids:
            df.loc[model_id, rxn_or_met] = 1

    dfs = [df]
    factors = ('organism', 'organism_id', 'template', 'method')
    labels = (organisms, organisms_id, templates, methods)

    for factor, label in zip(factors, labels):
        df = pd.DataFrame(data=label,
                          index=index,
                          columns=[factor])

        dfs.append(df)

    df = pd.concat(dfs, axis=1)

    return df


def cog_dataframe(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep='\t')
    df.index = df['organism']

    return df


def scaling(dataframe: pd.DataFrame,
            factors: Union[List[str], Tuple[str]],
            standard: bool = True,
            variance: bool = True) -> pd.DataFrame:
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
                 components: int = 2):
    x_mask = dataframe.columns[~dataframe.columns.isin(factors)]
    x = dataframe.loc[:, x_mask]
    y = dataframe.loc[:, factors]

    pca = PCA(n_components=components)

    pc = pca.fit_transform(x)

    columns = [f'PC {i + 1}' for i in range(components)]

    df = pd.DataFrame(data=pc, index=dataframe.index, columns=columns)

    return pd.concat([df, y], axis=1)


def plot_pca(workdir: str,
             dataframe: pd.DataFrame,
             pc1: str,
             pc2: str,
             factors: Union[List[str], Tuple[str]],
             content: str):

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    for factor in factors:

        fig = plt.figure(figsize=(8, 8))

        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(pc1, fontsize=15)
        ax.set_ylabel(pc2, fontsize=15)
        ax.set_title(f'{content} PCA', fontsize=20)

        labels = set(dataframe.loc[:, factor])

        colors = list(mcolors.TABLEAU_COLORS.keys())
        diff = len(labels) - len(colors)

        if diff > 0:
            colors += ['#1f6357', '#017374', '#0cb577', '#ff0789', '#afa88b']

        for label, color in zip(labels, colors):
            mask = dataframe.loc[:, factor] == label

            pc1_values = dataframe.loc[mask, pc1]
            pc2_values = dataframe.loc[mask, pc2]
            ax.scatter(pc1_values,
                       pc2_values,
                       c=color,
                       s=50)

            organisms_id = dataframe.loc[mask, 'organism_id']

            for pc1_pt, pc2_pt, annotation in zip(pc1_values, pc2_values, organisms_id):

                ax.annotate(annotation, (pc1_pt+1, pc2_pt-1))

        legend = ax.legend(labels, loc=(1.04, 0))
        ax.grid()
        file_name = f'{content}_{factor.title()}_{pc1}_{pc2}.png'.replace(' ', '_')
        file_path = os.path.join(workdir, file_name)
        fig.savefig(fname=file_path, bbox_extra_artists=(legend,), bbox_inches='tight')


def rxns_pca(models_dir: str, analysis_dir: str, filter_exchanges: bool):

    factors = ('organism', 'organism_id', 'template', 'method')

    analysis_models = read_models(models_dir)
    df = models_dataframe(analysis_models, reactions=True, filter_exchanges=filter_exchanges)
    df = scaling(df, factors=factors)
    pca = pca_analysis(df, factors=factors)
    plot_pca(workdir=analysis_dir, dataframe=pca, pc1='PC 1', pc2='PC 2',
             factors=('organism', 'template', 'method'), content='Reactions')


def mets_pca(models_dir: str, analysis_dir: str, filter_exchanges: bool):
    factors = ('organism', 'organism_id', 'template', 'method')

    analysis_models = read_models(models_dir)
    df = models_dataframe(analysis_models, reactions=False, filter_exchanges=filter_exchanges)
    df = scaling(df, factors=factors)
    pca = pca_analysis(df, factors=factors)
    plot_pca(workdir=analysis_dir, dataframe=pca, pc1='PC 1', pc2='PC 2',
             factors=('organism', 'template', 'method'), content='Metabolites')


def cog_pca(cog_analysis_file:str, analysis_dir: str):

    factors = ('domain', 'phylum')

    df = cog_dataframe(cog_analysis_file)
    df = scaling(df, factors)
    pca = pca_analysis(df, factors=factors, components=2)
    plot_pca(workdir=analysis_dir, dataframe=pca, pc1='PC 1', pc2='PC 2', factors=factors, content='COG')


if __name__ == '__main__':
    rxns_pca(os.path.join(os.getcwd(), 'models'),
             os.path.join(os.getcwd(), 'model_content_analysis'),
             filter_exchanges=True)
    mets_pca(os.path.join(os.getcwd(), 'models'),
             os.path.join(os.getcwd(), 'model_content_analysis'),
             filter_exchanges=True)
