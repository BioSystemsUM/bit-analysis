import os
import random
from collections import defaultdict

import pandas as pd

from comparative_analysis import read_bigg_models, get_organisms

random.seed(2)


def get_random_models(*species, organisms):
    randoms = {}

    for specie in species:

        randoms[specie] = []

        for _ in range(5):

            random_organisms = random.sample(organisms, 3)

            for random_organism in random_organisms:
                random_model = random.choice([model for model in random_organism.bigg_ids])

                randoms[specie].append((random_organism.name, random_model))

    return randoms


def write_random_report(workdir, models):
    dict_df = defaultdict(list)

    for organism, iterations in models.items():

        for organism_mame, bigg_model in iterations:
            dict_df[f'{organism}_name'].append(organism_mame)
            dict_df[f'{organism}_model'].append(bigg_model)

    df = pd.DataFrame.from_dict(dict_df)

    file_path = os.path.join(workdir, 'random_report.xlsx')
    df.to_excel(file_path)

    return df


if __name__ == '__main__':
    directory = os.path.join(os.getcwd(), 'comparative_func_analysis')

    bigg_models = read_bigg_models(workdir=directory)

    bigg_organisms = get_organisms(bigg_models, verbose=True)

    random_models = get_random_models('sth', 'mtb', 'xfa', organisms=bigg_organisms)

    final_df = write_random_report(directory, random_models)
    print()
