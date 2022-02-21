import re
from abc import ABC, abstractmethod

import pandas as pd


class XRefsConverter(ABC):

    def __init__(self, file_path):
        self.source_to_external_database_map = {}
        self.external_database_to_source_map = {}
        self.parse_xrefs_file(file_path)

    @property
    def source_to_external_database_map(self):
        return self._source_to_external_database_map

    @source_to_external_database_map.setter
    def source_to_external_database_map(self, value):
        self._source_to_external_database_map = value

    @property
    def external_database_to_source_map(self):
        return self._external_database_to_source_map

    @external_database_to_source_map.setter
    def external_database_to_source_map(self, value):
        self._external_database_to_source_map = value

    def convert(self, id_1, db_2):

        if id_1 in self.external_database_to_source_map:
            internal_id = self.external_database_to_source_map[id_1]

            if db_2 in self.source_to_external_database_map[internal_id]:
                return self.source_to_external_database_map[internal_id][db_2]

        return None

    @abstractmethod
    def parse_xrefs_file(self, file_path):
        raise NotImplementedError


class MetaNetXCompoundsConverter(XRefsConverter):

    def parse_xrefs_file(self, file_path):
        df = pd.read_csv(file_path)

        for i, row in df.iterrows():
            external_id = row["#source"]
            internal_id = row["ID"]
            if re.match("kegg\.compound:.*", external_id):
                external_id = external_id.replace("kegg.compound:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"kegg": [external_id]}

                elif "kegg" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["kegg"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["kegg"].append(external_id)


            elif re.match("bigg\.metabolite:.*", external_id):
                external_id = external_id.replace("bigg.metabolite:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"bigg": [external_id]}

                elif "bigg" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["bigg"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["bigg"].append(external_id)

            elif re.match("metacyc\.compound:.*", external_id):
                external_id = external_id.replace("metacyc.compound:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"metacyc": [external_id]}

                elif "metacyc" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["metacyc"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["metacyc"].append(external_id)

            elif re.match("seed\.compound:.*", external_id):
                external_id = external_id.replace("seed.compound:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"seed": [external_id]}

                elif "seed" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["seed"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["seed"].append(external_id)


class MetaNetXReactionsConverter(XRefsConverter):

    def parse_xrefs_file(self, file_path):
        df = pd.read_csv(file_path)

        for i, row in df.iterrows():
            external_id = row["#source"]
            internal_id = row["ID"]
            if re.match("kegg\.reaction:.*", external_id):
                external_id = external_id.replace("kegg.reaction:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"kegg": [external_id]}

                elif "kegg" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["kegg"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["kegg"].append(external_id)

            elif re.match("bigg\.reaction:.*", external_id):
                external_id = external_id.replace("bigg.reaction:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"bigg": [external_id]}

                elif "bigg" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["bigg"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["bigg"].append(external_id)

            elif re.match("metacyc\.reaction:.*", external_id):
                external_id = external_id.replace("metacyc.reaction:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"metacyc": [external_id]}

                elif "metacyc" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["metacyc"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["metacyc"].append(external_id)

            elif re.match("seed\.reaction:.*", external_id):
                external_id = external_id.replace("seed.reaction:", "")
                self.external_database_to_source_map[external_id] = internal_id

                if internal_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[internal_id] = {"seed": [external_id]}

                elif "seed" not in self.source_to_external_database_map[internal_id]:
                    self.source_to_external_database_map[internal_id]["seed"] = [external_id]

                else:
                    self.source_to_external_database_map[internal_id]["seed"].append(external_id)

class SeedReactionsConverter(XRefsConverter):

    def parse_xrefs_file(self, file_path):
        df = pd.read_csv(file_path)

        for i, row in df.iterrows():
            external_id = row["External ID"]
            seed_id = row["ModelSEED ID"]
            database = row["Source"]
            if database == "KEGG":
                self.external_database_to_source_map[external_id] = seed_id

                if seed_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[seed_id] = {"kegg": [external_id]}

                elif "kegg" not in self.source_to_external_database_map[seed_id]:
                    self.source_to_external_database_map[seed_id]["kegg"] = [external_id]

                else:
                    self.source_to_external_database_map[seed_id]["kegg"].append(external_id)

            elif database == "BiGG":
                self.external_database_to_source_map[external_id] = seed_id

                if seed_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[seed_id] = {"bigg": [external_id]}

                elif "bigg" not in self.source_to_external_database_map[seed_id]:
                    self.source_to_external_database_map[seed_id]["bigg"] = [external_id]

                else:
                    self.source_to_external_database_map[seed_id]["bigg"].append(external_id)

            elif database == "MetaCyc":
                self.external_database_to_source_map[external_id] = seed_id

                if seed_id not in self.source_to_external_database_map:
                    self.source_to_external_database_map[seed_id] = {"metacyc": [external_id]}

                elif "metacyc" not in self.source_to_external_database_map[seed_id]:
                    self.source_to_external_database_map[seed_id]["metacyc"] = [external_id]

                else:
                    self.source_to_external_database_map[seed_id]["metacyc"].append(external_id)

            # ModelSEED IDs

            if seed_id not in self.external_database_to_source_map:
                self.external_database_to_source_map[seed_id] = seed_id

            if seed_id not in self.source_to_external_database_map:
                self.source_to_external_database_map[seed_id] = {"seed": [seed_id]}

            elif "seed" not in self.source_to_external_database_map[seed_id]:
                self.source_to_external_database_map[seed_id]["seed"] = [seed_id]

            else:
                if seed_id not in self.source_to_external_database_map[seed_id]["seed"]:
                    self.source_to_external_database_map[seed_id]["seed"].append(seed_id)
