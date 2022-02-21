import json
import re
from abc import ABC, abstractmethod
from utils import ReconstructionTool

# TODO : identify automatically the database version
import cobra.io


def read_sbml_into_cobra_model(file_path, database_version, reconstruction_tool):
    new_model = CobraModel()

    new_model.model = cobra.io.read_sbml_model(file_path)
    new_model.database_version = database_version
    new_model.reconstruction_tool = reconstruction_tool

    return new_model


class Report(ABC):

    @abstractmethod
    def save_to_json(self, file_path):
        raise NotImplementedError


class ModelComparisonReport(Report):

    def save_to_json(self, file_path):
        pass


class FormatConversionReport(Report):

    def __init__(self):
        self._non_convertable = None
        self._convertable = None
        self._conversion_map = {}

    @property
    def conversion_map(self):
        return self._conversion_map

    @conversion_map.setter
    def conversion_map(self, value):
        self.conversion_map = value

    @property
    def convertable(self):
        return self._convertable

    @convertable.setter
    def convertable(self, value: dict):
        if isinstance(value, dict):
            self._convertable = value
            self.conversion_map["convertable"] = self.convertable

        else:
            raise TypeError

    @property
    def non_convertable(self):
        return self._non_convertable

    @non_convertable.setter
    def non_convertable(self, value):
        if isinstance(value, list):
            self._non_convertable = value
            self.conversion_map["non_convertable"] = self.non_convertable

        else:
            raise TypeError

    def save_to_json(self, file_path):
        json.dump(self.conversion_map, file_path)


class Model(ABC):

    def __init__(self):
        self._reaction_converter = None
        self._model = None
        self._database_version = None
        self._metabolite_converter = None
        self._reconstruction_tool = None

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        self._model = value

    @property
    def database_format(self):
        return self._database_version

    @database_format.setter
    def database_format(self, value):
        self._database_format = value

    @property
    def reconstruction_tool(self):
        return self._reconstruction_tool

    @reconstruction_tool.setter
    def reconstruction_tool(self, value):
        self._reconstruction_tool = value

    # TODO : put database into enumerators
    @abstractmethod
    def get_metabolites_other_version(self, database):
        raise NotImplementedError

    @abstractmethod
    def get_reactions_other_version(self, database):
        raise NotImplementedError

    @abstractmethod
    def convert_model_into_other_database(self, database):
        raise NotImplementedError

    @property
    def metabolite_converter(self):
        return self._metabolite_converter

    @metabolite_converter.setter
    def metabolite_converter(self, value):
        self._metabolite_converter = value

    @property
    def reaction_converter(self):
        return self._reaction_converter

    @reaction_converter.setter
    def reaction_converter(self, value):
        self._reaction_converter = value


class CobraModel(Model):

    def convert_model_into_other_database(self, database):
        pass

    def get_metabolites_other_version(self, database) -> FormatConversionReport:
        convertable = {}
        not_convertable = []

        regex = re.compile("_[a-zA-Z]+$")

        for metabolite in self.model.metabolites:
            metabolite_id = metabolite.id

            # if self.reconstruction_tool == ReconstructionTool.MERLIN:
            #     metabolite_id = str(metabolite.id).split("__")[0]
            # elif self.reconstruction_tool == ReconstructionTool.MODELSEED or \
            #         self.reconstruction_tool == ReconstructionTool.CARVEME or \
            #         self.reconstruction_tool == ReconstructionTool.AUREME or \
            #         self.reconstruction_tool == ReconstructionTool.PATHWAYTOOLS:
            #     metabolite_id = str(metabolite.id).split("_")[0]

            metabolite_id = re.sub(regex, "", metabolite_id)
            metabolite_list = self.metabolite_converter.convert(metabolite_id, database)
            if metabolite_list:
                convertable[metabolite_id] = metabolite_list
            else:
                not_convertable.append(metabolite_id)

        report = FormatConversionReport()
        report.convertable = convertable
        report.non_convertable = not_convertable

        return report

    def get_reactions_other_version(self, database):
        convertable = {}
        not_convertable = []

        for reaction in self.model.reactions:

            reaction_id = reaction.id

            if self.reconstruction_tool == ReconstructionTool.MERLIN:
                reaction_id = str(reaction.id).split("__")[0]
            elif self.reconstruction_tool == ReconstructionTool.MODELSEED:
                reaction_id = str(reaction.id).split("_")[0]

            reaction_list = self.reaction_converter.convert(reaction_id, database)
            if reaction_list:
                convertable[reaction_id] = reaction_list
            else:
                not_convertable.append(reaction_id)

        report = FormatConversionReport()
        report.convertable = convertable
        report.non_convertable = not_convertable

        return report
