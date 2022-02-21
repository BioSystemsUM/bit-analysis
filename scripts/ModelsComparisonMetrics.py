import cobra
import core
from utils import Utils, Type

import matplotlib
matplotlib.use('Agg')

import seaborn
seaborn.set_style("darkgrid")

import venn
import pandas

class ModelsComparisonMetrics:

    def __init__(self, models):
        self.models = models
        self.metabolite_sets = {}
        self.reaction_sets = {}
        self.gene_sets = {}

    def set_models(self, models) -> None:
        self.models = models

    def set_metabolite_sets(self, metabolite_sets) -> None:
        self.metabolite_sets = metabolite_sets

    def set_reaction_sets(self, reaction_sets) -> None:
        self.reaction_sets = reaction_sets

    def set_gene_sets(self, gene_sets) -> None:
        self.gene_sets = gene_sets

    def __gene_sets(self) -> None:
        for model_name in self.models.keys():
            gene_set = []
            model = self.models[model_name]
            for gene in model.genes:
                gene_set.append(gene.id)
            self.gene_sets[model_name] = set(gene_set)

    def __reaction_sets(self) -> None:
        for model_name in self.models.keys():
            reaction_set = []
            model = self.models[model_name]
            for reaction in model.reactions:
                reaction_set.append(reaction.id)
            self.reaction_sets[model_name] = set(reaction_set)

    def __metabolites_sets(self) -> None:
        for model_name in self.models.keys():
            metabolites_set = []
            model = self.models[model_name]
            for metabolite in model.metabolites:
                metabolites_set.append(metabolite.id)
            self.metabolite_sets[model_name] = set(metabolites_set)

    def calculate_jaccard_distances(self, reference_model_name, type) -> dict:

        jaccard_distances = {}

        if reference_model_name not in self.models.keys():
            raise Exception("There is no model called " + str(reference_model_name))
        else:

            if type == Type.GENES:

                if self.gene_sets == {}:
                    self.__gene_sets()

                for model_name in self.gene_sets.keys():
                    if model_name != reference_model_name:
                        jaccard_distance = Utils.jaccard_distance(self.gene_sets[model_name],
                                                                  self.gene_sets[reference_model_name])
                        jaccard_distances[model_name] = jaccard_distance

            elif type == Type.REACTIONS:

                if self.reaction_sets == {}:
                    self.__reaction_sets()

                for model_name in self.reaction_sets.keys():
                    if model_name != reference_model_name:
                        jaccard_distance = Utils.jaccard_distance(self.reaction_sets[model_name],
                                                                  self.reaction_sets[reference_model_name])
                        jaccard_distances[model_name] = jaccard_distance

                        disjunction= self.reaction_sets[reference_model_name] - self.reaction_sets[model_name]

                        file =open('results/' + model_name + "_" + self.models['carveme model'].id +"_disjunction.txt", "w")
                        for reaction in disjunction:
                            file.write(reaction + "\n")
                        file.close()

            elif type == Type.METABOLITES:

                if self.metabolite_sets == {}:
                    self.__metabolites_sets()

                for model_name in self.metabolite_sets.keys():
                    if model_name != reference_model_name:
                        jaccard_distance = Utils.jaccard_distance(self.metabolite_sets[model_name],
                                                                  self.metabolite_sets[reference_model_name])
                        jaccard_distances[model_name] = jaccard_distance

        return jaccard_distances

    def calculate_ratios(self, reference_model_name, type) -> dict:

        ratios = {}

        if reference_model_name not in self.models.keys():
            raise Exception("There is no model called " + str(reference_model_name))

        else:

            if type == Type.GENES:

                if self.gene_sets == {}:
                    self.__gene_sets()

                for model_name in self.gene_sets.keys():
                    if model_name != reference_model_name:
                        ratio = Utils.ratio(self.gene_sets[model_name], self.gene_sets[reference_model_name])
                        ratios[model_name] = ratio

            elif type == Type.REACTIONS:

                if self.reaction_sets == {}:
                    self.__reaction_sets()

                for model_name in self.reaction_sets.keys():
                    if model_name != reference_model_name:
                        ratio = Utils.ratio(self.reaction_sets[model_name], self.reaction_sets[reference_model_name])
                        ratios[model_name] = ratio

            elif type == Type.METABOLITES:

                if self.metabolite_sets == {}:
                    self.__metabolites_sets()

                for model_name in self.metabolite_sets.keys():
                    if model_name != reference_model_name:
                        ratio = Utils.ratio(self.metabolite_sets[model_name],
                                            self.metabolite_sets[reference_model_name])
                        ratios[model_name] = ratio

        return ratios

    @staticmethod
    def draw_dot_plot(jaccard_distances, ratios, output_path):

        listofzeros = [0.00000] * len(jaccard_distances.keys())

        df = pandas.DataFrame({
            "model": jaccard_distances.keys(),
            "Jaccard distance": listofzeros,
            "Ratio": listofzeros
        })

        for index, row in df.iterrows():
            model = row["model"]
            df._set_value(index, "Jaccard distance", jaccard_distances[model])
            df._set_value(index, "Ratio", ratios[model])


        matplotlib.pyplot.clf()
        # fig = matplotlib.pyplot.figure(figsize=(25, 10))

        #seaborn.cubehelix_palette(start=.5, rot=-.75, as_cmap=True)

        # matplotlib.pyplot.subplot(121)
        # dotplot0 = seaborn.scatterplot(data=df, x="jaccard distance", y="ratio", hue="model", size="model", palette="viridis", legend=False)
        # dotplot0.set_xlim([0, 1])
        # dotplot0.set_ylim([0, 1])
        #
        # matplotlib.pyplot.subplot(122)
        dotplot = seaborn.scatterplot(data=df, x="Jaccard distance", y="Ratio", hue="model", style="model", size="model",
                                      sizes=(50, 100), palette="crest", legend="full")
        dotplot.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), facecolor="inherit", edgecolor="inherit", ncol=2,
                       fontsize='small')

        # xticks = [0.6500]
        # while (xticks[-1] < 0.7800):
        #     xticks.append(xticks[-1]+0.05)
        #
        # yticks = [0.4000]
        # while (yticks[-1] < 0.7000):
        #     yticks.append(yticks[-1]+0.05)
        #
        # dotplot.set_xlim([0.6500, 0.7800])
        # dotplot.set_xticks(xticks)
        #
        # dotplot.set_ylim([0.4000, 0.7000])
        # dotplot.set_yticks(yticks)

        dotplot.tick_params(axis='x', rotation=20, labelsize="x-small")
        dotplot.tick_params(axis='y', rotation=20, labelsize="x-small")

        #figure = fig.get_figure()
        figure = dotplot.get_figure()
        figure.savefig(output_path, dpi=1200, bbox_inches="tight")

    def draw_venn_diagram(self, type, output_path) -> None:

        data = set()

        if type == Type.GENES:

            if len(self.gene_sets.keys()) > 6:
                raise Exception("Venn diagram construction only accepts a maximum of 6 models")

            if self.gene_sets == {}:
                self.__gene_sets()

            data = self.gene_sets

        elif type == Type.REACTIONS:

            if len(self.reaction_sets.keys()) > 6:
                raise Exception("Venn diagram construction only accepts a maximum of 6 models")

            if self.reaction_sets == {}:
                self.__reaction_sets()

            data = self.reaction_sets

        elif type == Type.METABOLITES:

            if len(self.metabolite_sets.keys()) > 6:
                raise Exception("Venn diagram construction only accepts a maximum of 6 models")

            if self.metabolite_sets == {}:
                self.__metabolites_sets()

            data = self.metabolite_sets

        matplotlib.pyplot.clf()

        if len(data.keys()) == 6:
            venndiagram = venn.pseudovenn(data, fmt="{percentage:.1f}%", cmap="crest_r", fontsize="small",
                                          legend_loc="upper left", hint_hidden=False)
        else:
            venndiagram = venn.venn(data, fmt="{percentage:.1f}%", cmap="crest_r", fontsize="small",
                                    legend_loc="upper left", hint_hidden=False)

        venndiagram.legend(data.keys(), loc='lower center', bbox_to_anchor=(0.5, -0.10), ncol=2, fontsize='small',
                           facecolor="inherit", edgecolor="inherit")
        figure = venndiagram.get_figure()
        figure.savefig(output_path, dpi=1200, bbox_inches="tight")

    @staticmethod
    def draw_barplot(execution_time_values, output_path):

        seaborn.set_context(rc={'patch.linewidth': 0.0})

        barplot = seaborn.barplot(data=execution_time_values, x="Duration (seconds)", y="Organism", hue="Tool / Sensitive mode",
                                  palette="crest")
        barplot.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), ncol=1, fontsize='small', facecolor="inherit", edgecolor="inherit")
        barplot.set_ylabel("")

        xticks = []

        for i in range(0, 18000, 1800):
            xticks.append(i)

        barplot.set_xticks(xticks)
        barplot.set_xlim([0, 18000])

        barplot.set_yticklabels(
            ["Bordetella\nPertussis", "Komagataeibacter\nSucrofermentans", "Lactobacillus\nPlantarum",
             "Pseudomonas\nPutida"])

        barplot.tick_params(axis='x', rotation=20, labelsize="x-small")
        barplot.tick_params(axis='y', rotation=20, labelsize="x-small")

        figure = barplot.get_figure()
        figure.savefig(output_path, dpi=1200, bbox_inches="tight")

    def write_reaction_genes_metabolites_into_csv(self, organism, output_path_non_curated_models, output_path_curated_models) -> None:

        organisms_list = []
        curated_organisms_list = []
        tools_list = []
        genes_count = []
        curated_genes_count = []
        reactions_count = []
        curated_reactions_count = []
        metabolites_count = []
        curated_metabolites_count = []
        data = {}
        curated_data = {}

        for model_name in self.models.keys():

            model = self.models[model_name]

            organism = organism
            tool = ""

            geneSet = self.gene_sets[model_name]
            reactionSet = self.reaction_sets[model_name]
            metabolitesSet = self.metabolite_sets[model_name]

            if isinstance(model, cobra.core.model.Model):

                if "carveme" in str(model_name):
                    tool = "carveme"
                elif "all" in str(model_name):
                    tool = "BIT"
                elif "selected" in str(model_name):
                    tool = "BIT"
                elif "random" in str(model_name):
                    tool = "BIT"
                elif "reference" in str(model_name):
                    tool = "reference"
                elif "curated" in str(model_name):
                    tool = "reference"

                organisms_list.append(organism)
                tools_list.append(tool)
                genes_count.append(len(geneSet))
                reactions_count.append(len(reactionSet))
                metabolites_count.append(len(metabolitesSet))

                data = {"Organism": organisms_list,
                        "Tool": tools_list,
                        "Number of Genes": genes_count,
                        "Number of Reactions": reactions_count,
                        "Number of Metabolites": metabolites_count
                        }

            elif isinstance(model, core.Model):

                curated_organisms_list.append(organism)

                curated_genes_count.append(len(geneSet))
                curated_reactions_count.append(len(reactionSet))
                curated_metabolites_count.append(len(metabolitesSet))

                curated_data = {"Organism": curated_organisms_list,
                                "Number of Genes": curated_genes_count,
                                "Number of Reactions": curated_reactions_count,
                                "Number of Metabolites": curated_metabolites_count
                                }

        if data != {}:
            dataframe = pandas.DataFrame(data)
            dataframe.to_csv(output_path_non_curated_models, sep="\t", index=False)

        if curated_data != {}:
            curated_dataframe = pandas.DataFrame(curated_data)
            curated_dataframe.to_csv(output_path_curated_models, sep="\t", index=False)

    @staticmethod
    def write_ratios_and_jaccard_distances(organism, jaccard_distances, ratios, output_path):

        organisms_list = []
        tools_list = []
        jaccard_distances_list = []
        ratios_list = []

        for model_name in jaccard_distances.keys():

            tool = ""

            if "carveme" in str(model_name):
                tool = "carveme"
            elif "blast" in str(model_name):
                tool = "blast"
            elif "diamond fast" in str(model_name):
                tool = "diamond fast"
            elif "diamond sensitive" in str(model_name):
                tool = "diamond sensitive"
            elif "diamond very-sensitive" in str(model_name):
                tool = "diamond very-sensitive"
            elif "diamond ultra-sensitive" in str(model_name):
                tool = "diamond ultra-sensitive"

            organisms_list.append(organism)
            tools_list.append(tool)
            jaccard_distances_list.append(jaccard_distances[model_name])
            ratios_list.append(ratios[model_name])

        jaccard_distances_list_round = [round(elem, 4) for elem in jaccard_distances_list]
        ratios_list_round = [round(elem, 4) for elem in ratios_list]

        data = {"Organism": organisms_list,
                "Tool": tools_list,
                "Jaccard distance": jaccard_distances_list_round,
                "Ratio": ratios_list_round
                }

        dataframe = pandas.DataFrame(data)
        dataframe.to_csv(output_path, sep="\t", index=False)

