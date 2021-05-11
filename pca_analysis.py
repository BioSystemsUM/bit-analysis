import os
from collections import defaultdict
from typing import List

import pandas as pd
from cobra import Model
from cobra.io import read_sbml_model
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
            cobra_model = read_sbml_model(model_path)
            cobra_models.append(cobra_model)

    return cobra_models


def features_dataframe(models: List[Model], model_feature: str = 'reactions') -> pd.DataFrame:

    features_lookup = defaultdict(list)

    index = []
    for model in models:

        index.append(model.id)

        model_attr = getattr(model, model_feature.lower(), [])
        for feature in model_attr:

            features_lookup[feature.id].append(model.id)

    data = [[0] * len(features_lookup)] * len(index)
    df = pd.DataFrame(data=data,
                      index=index,
                      columns=features_lookup.keys())

    for rxn, models_ids in features_lookup.items():

        for model_id in models_ids:
            df.loc[model_id, rxn] = 1

    return df


def scaling(dataframe: pd.DataFrame, standard=True, variance=True):

    out_df = dataframe

    if standard:
        t_df = out_df.T

        scalier = StandardScaler(with_mean=False, with_std=False)
        scaled = scalier.fit_transform(t_df)

        out_df = pd.DataFrame(scaled, columns=dataframe.index, index=dataframe.columns)
        out_df = out_df.T

    if variance:
        scalier = VarianceThreshold()
        scaled = scalier.fit(out_df)
        out_df = out_df.iloc[:, scaled.get_support(indices=True)]

    return out_df


def pca(dataframe: pd.DataFrame):

    pass


def annotate():
    pass


if __name__ == '__main__':

    directory = os.path.join(os.getcwd(), 'models')
    analysis_models = read_models(directory)
    rxns_df = features_dataframe(analysis_models)
    scaled_df = scaling(rxns_df)
    print()