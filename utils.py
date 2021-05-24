import os
import re
from collections import namedtuple
from typing import List, Union

from matplotlib import colors as mcolors

from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model

M_COLORS = list(mcolors.TABLEAU_COLORS.keys())
COLORS = M_COLORS + ['#1f6357', '#017374', '#0cb577', '#ff0789', '#afa88b']

MARKERS = ['o', '^', '*', '8', 's', 'p']

metabolite_pattern = re.compile(r'(__[A-Za-z](?!.))')
extracellular_metabolite_pattern = re.compile(r'(__[A-Za-z]_b(?!.))')
reaction_pattern = re.compile(r'(__[A-Za-z]+(?!.))')

ModelAnalysis = namedtuple('ModelAnalysis',
                           field_names=('model',
                                        'model_id',
                                        'organism',
                                        'organism_id',
                                        'template',
                                        'method'),
                           defaults=(Model(),
                                     'Mtuberculosis_all_permissive'
                                     'Mycobacterium tuberculosis',
                                     'Mtub'
                                     'all',
                                     'permissive'))


def parse_organism_annotation(model_annotation):
    if len(model_annotation) > 1:

        organism = model_annotation[0]

        if 'tuberculosis' in organism.lower():
            return 'Mycobacterium tuberculosis'

        elif 'thermophilus' in organism.lower():
            return 'Streptococcus thermophilus'

        elif 'fastidiosa' in organism.lower():
            return 'Xylella fastidiosa'

    return


def parse_organism_id(organism_name):
    genus_species = organism_name.split()
    genus = genus_species[0]
    species = genus_species[1]

    return f'{genus[0]}{species[0:3]}'


def parse_template_annotation(model_annotation):
    if len(model_annotation) > 1:

        template = model_annotation[1]

        if 'all' in template.lower():
            return 'all', None

        elif 'random' in template.lower():
            return 'random', template.lower().replace('random', '')

        elif 'select' in template.lower():
            return 'select', None

        elif 'carveme' in template.lower():
            return 'carveme', None

    return


def parse_method_annotation(model_annotation):
    if len(model_annotation) > 1:

        method = model_annotation[2]

        if 'permissive' in method:
            return 'permissive'

        elif 'restrictive' in method:
            return 'restrictive'

        elif 'carveme' in method.lower():
            return 'carveme'

    return


def read_models(workdir: str) -> List[ModelAnalysis]:
    if not os.path.exists(workdir):
        raise OSError(f'{workdir} does not exist')

    models = os.listdir(workdir)

    models_analysis = []

    for model_name in models:

        if model_name.endswith('.xml') or model_name.endswith('.sbml'):
            model_path = os.path.join(workdir, model_name)
            model: Model = read_sbml_model(model_path)

            model_id = model_name.replace('model_', '').replace('.xml', '').replace('.sbml', '')

            model_annotation = model_id.split('_')

            organism = parse_organism_annotation(model_annotation)
            organism_id = parse_organism_id(organism)
            template, random_id = parse_template_annotation(model_annotation)

            if template == 'random':
                organism_id = f'{organism_id}{random_id}'

            method = parse_method_annotation(model_annotation)

            model.id = model_id

            model_analysis = ModelAnalysis(model=model,
                                           model_id=model_id,
                                           organism=organism,
                                           organism_id=organism_id,
                                           template=template,
                                           method=method)

            models_analysis.append(model_analysis)

    return models_analysis


def parse_metabolite_id(metabolite_id: str) -> str:
    """
    merlin regularly appends either __e, __e_b, ... to a BiGG metabolite id.
    This function clips such suffixes using regex.
    :param metabolite_id: str, string for the metabolite
    :return: str, the real BiGG metabolite identifier
    """

    if re.search(metabolite_pattern, metabolite_id):
        return metabolite_pattern.sub('', metabolite_id)

    elif re.search(extracellular_metabolite_pattern, metabolite_id):
        return extracellular_metabolite_pattern.sub('_b', metabolite_id)

    else:
        return metabolite_id


def parse_reaction_id(reaction_id: str) -> str:
    """
    merlin regularly appends either __cytop, ... to a BiGG reaction id.
    This function clips such suffixes using regex.
    :param reaction_id: str, string for the metabolite
    :return: str, the real BiGG metabolite identifier
    """

    if re.search(reaction_pattern, reaction_id):
        return reaction_pattern.sub('', reaction_id)

    else:
        return reaction_id


def parse_reaction(reaction: Reaction, filter_boundaries: bool = True) -> Union[str, None]:
    """
    Parsing reaction to the correct reaction BiGG identifier.
    It returns None if the reaction is boundary and filter_boundaries is True
    :param reaction: Reaction, a cobra reaction object
    :param filter_boundaries: bool, Whether boundary reactions should be filtered
    :return: str or None, the correct reaction BiGG identifier
    """

    if filter_boundaries:

        is_exchange = False

        if reaction.boundary:
            is_exchange = True

        else:

            for met in reaction.metabolites:
                if met.id.endswith('_b'):
                    is_exchange = True
                    break

        if is_exchange:
            return

    return parse_reaction_id(reaction.id)


def parse_metabolite(metabolite: Metabolite, filter_boundaries: bool = True) -> Union[str, None]:
    """
    Parsing metabolite to the correct metabolite BiGG identifier.
    It returns None if the metabolite is boundary and filter_boundaries is True
    :param metabolite: Metabolite, a cobra metabolite object
    :param filter_boundaries: bool, Whether boundary metabolites should be filtered
    :return: str or None, the correct metabolite BiGG identifier
    """

    if filter_boundaries and metabolite.id.endswith('_b'):
        return

    return parse_metabolite_id(metabolite.id)


def get_explained_variance_idx(component):
    if '1' in component:
        return 0

    elif '2' in component:
        return 1

    elif '3' in component:
        return 2

    return 3