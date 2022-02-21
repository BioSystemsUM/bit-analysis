import re
from enum import Enum
import cobra.io.sbml

class Utils:

    @staticmethod
    def get_metabolite_ids(model_path):
        model = cobra.io.sbml.read_sbml_model(model_path)
        metabolite_list = []

        for metabolite in model.metabolites:
            metabolite_id = re.sub("_[a-z]*$", "", metabolite.id)
            metabolite_list.append(metabolite_id)

        metabolite_set = set(metabolite_list)

        with open("metabolites_list_%s.txt" % model.name, "w") as file:
            for met_id in metabolite_set:
                file.write(met_id + "\n")

    @staticmethod
    def jaccard_distance(generated_elements, reference_elements):

        jd = 0

        if type(generated_elements) != set:
            raise TypeError("generated_elements should be a set instead of " + str(type(generated_elements)))

        elif type(reference_elements) != set:
            raise TypeError("reference_elements should be a set instead of " + str(type(reference_elements)))

        else:
            intersection = generated_elements & reference_elements
            union = generated_elements | reference_elements
            jd = 1 - (abs(len(intersection))/abs(len(union)))

        return jd

    @staticmethod
    def ratio(generated_elements, reference_elements):

        ratio = 0

        if type(generated_elements) != set:
            raise TypeError("generated_elements should be a set instead of " + str(type(generated_elements)))

        elif type(reference_elements) != set:
            raise TypeError("reference_elements should be a set instead of " + str(type(reference_elements)))

        else:

            intersection = generated_elements & reference_elements
            difference = generated_elements - reference_elements
            ratio = abs(len(intersection)) / abs(len(difference))

        return ratio


class Type(Enum):
    GENES = 1
    REACTIONS = 2
    METABOLITES = 3

class ReconstructionTool(Enum):
    AUTOKEGGREC = 1
    AUREME = 2
    CARVEME = 3
    MERLIN = 4
    MODELSEED = 5
    PATHWAYTOOLS = 6
    RAVEN = 7


metabolite_pattern = re.compile(r'(__[A-Za-z](?!.))')
extracellular_metabolite_pattern = re.compile(r'(__[A-Za-z]_b(?!.))')
reaction_pattern = re.compile(r'(__[A-Za-z]+(?!.))')

def parse_reaction_id(reaction):
    """
    merlin regularly appends either __cytop, ... to a BiGG reaction id.
    This function clips such suffixes using regex.
    :param reaction_id: str, string for the metabolite
    :return: str, the real BiGG metabolite identifier
    """

    if re.search(reaction_pattern, reaction.id):
        return reaction_pattern.sub('', reaction.id)

    else:
        return reaction


def parse_reaction(reaction, filter_boundaries: bool = True):
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

    return parse_reaction_id(reaction)


if __name__ == "__main__":
    Utils.get_metabolite_ids("PseudomonasPutida-blast-model.xml")
