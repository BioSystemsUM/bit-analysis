import matplotlib
import seaborn

from xrefs_converters import MetaNetXCompoundsConverter, MetaNetXReactionsConverter, SeedReactionsConverter
from ModelsComparisonMetrics import ModelsComparisonMetrics
from core import read_sbml_into_cobra_model
from utils import Utils, Type, ReconstructionTool
from cobra.io.sbml import read_sbml_model
import pandas
from Bio import SeqIO, GenBank


def getModels(organism) -> dict:

    # comparison models
    aureme_model = read_sbml_into_cobra_model(file_path="../Models/AuReMe/" + str(organism) + ".sbml",
                                              database_version="bigg",
                                              reconstruction_tool=ReconstructionTool.AUREME)

    autokeggrec_model = read_sbml_into_cobra_model(file_path="../Models/AutoKEGGRec/" + str(organism) + ".xml",
                                                   database_version="kegg",
                                                   reconstruction_tool=ReconstructionTool.AUTOKEGGREC)

    carveme_model = read_sbml_into_cobra_model(file_path="../Models/CarveMe/" + str(organism) + ".xml",
                                               database_version="bigg",
                                               reconstruction_tool=ReconstructionTool.CARVEME)

    modelseed_model = read_sbml_into_cobra_model(file_path="../Models/ModelSEED/" + str(organism) + ".sbml",
                                                 database_version="seed",
                                                 reconstruction_tool=ReconstructionTool.MODELSEED)

    pathwaytools_model = read_sbml_into_cobra_model(file_path="../Models/PathwayTools/" + str(organism) + ".xml",
                                                    database_version="metacyc",
                                                    reconstruction_tool=ReconstructionTool.PATHWAYTOOLS)

    raven_model = read_sbml_into_cobra_model(file_path="../Models/RAVEN/" + str(organism) + ".xml",
                                             database_version="kegg",
                                             reconstruction_tool=ReconstructionTool.RAVEN)

    merlin_model = read_sbml_into_cobra_model(file_path="../Models/Merlin/" + str(organism) + ".xml",
                                              database_version="kegg",
                                              reconstruction_tool=ReconstructionTool.MERLIN)

    # reference model
    reference_model = read_sbml_into_cobra_model(file_path="../Models/Manually_curated/" + str(organism) + ".xml",
                                                 database_version="bigg",
                                                 reconstruction_tool=None)


    models = {"aureme model": aureme_model,
              "autokeggrec model": autokeggrec_model,
              "carveme model": carveme_model,
              "modelseed model": modelseed_model,
              "pathwaytools model": pathwaytools_model,
              "raven model": raven_model,
              "merlin model": merlin_model,
              "reference model": reference_model
              }

    return models

def getReactionSets(models) -> dict:

    reaction_sets = {}

    for model_name in models.keys():
        reaction_set = []
        model = models[model_name]
        #reactions_converter = MetaNetXReactionsConverter("xrefs_files/reac_xrefs_reduced.csv")
        reactions_converter = SeedReactionsConverter("xrefs_files/reac_xrefs_reduced_from_seed_aliases.csv")
        model.reaction_converter = reactions_converter
        report = model.get_reactions_other_version("seed")
        convertable_reactions = report.convertable
        for convertable_reaction in convertable_reactions:
            kegg_reactions = convertable_reactions[convertable_reaction]

            if len(kegg_reactions) > 1:
                dic = {}
                for kegg_reaction in kegg_reactions:
                    dic[kegg_reaction] = 0
                    for model in reaction_sets.keys():
                        if kegg_reaction in reaction_sets[model]:
                            dic[kegg_reaction] = int(dic[kegg_reaction]) + 1

                dic = dict(sorted(dic.items(), key=lambda item: item[1]))
                reaction_set.append(list(dic.keys())[-1])
            else:
                for kegg_reaction in kegg_reactions:
                    if kegg_reactions not in reaction_set:
                        reaction_set.append(kegg_reaction)
        reaction_sets[model_name] = set(reaction_set)
    return reaction_sets

def getGeneSets(models) -> dict:

    gene_sets = {}

    for model_name in models.keys():
        gene_set = []
        model = models[model_name]
        for gene in model.model.genes:
            gene_set.append(gene.id)
        gene_sets[model_name] = set(gene_set)
    return gene_sets

def resultsGeneration (organism, genes_dotplot_path, reactions_dotplot_path) -> None:

    models = getModels(organism=organism)
    modelsComparisonMetrics = ModelsComparisonMetrics(models=models)

    gene_sets = getGeneSets(models=models)
    modelsComparisonMetrics.set_gene_sets(gene_sets=gene_sets)

    reaction_sets = getReactionSets(models=models)
    modelsComparisonMetrics.set_reaction_sets(reaction_sets=reaction_sets)

    for modelName in models.keys():

        model = models[modelName]
        geneSet = gene_sets[modelName]
        reactionSet = reaction_sets[modelName]

        print(modelName)
        print(len(model.model.genes))
        print(len(geneSet))
        print(len(model.model.reactions))
        print(len(reactionSet))
        print()

    # Genes

    jaccard_distances_genes = modelsComparisonMetrics.calculate_jaccard_distances(
        reference_model_name="reference model",
        type=Type.GENES)
    ratios_genes = modelsComparisonMetrics.calculate_ratios(
        reference_model_name="reference model",
        type=Type.GENES)

    gene_sets.pop("reference model")
    modelsComparisonMetrics.set_gene_sets(gene_sets=gene_sets)

    modelsComparisonMetrics.draw_dot_plot(jaccard_distances=jaccard_distances_genes,
                                          ratios=ratios_genes,
                                          output_path=genes_dotplot_path)

    # Reactions

    jaccard_distances_reactions = modelsComparisonMetrics.calculate_jaccard_distances(reference_model_name="reference model",
                                                                                      type=Type.REACTIONS)
    ratios_reactions = modelsComparisonMetrics.calculate_ratios(reference_model_name="reference model",
                                                                type=Type.REACTIONS)

    reaction_sets.pop("reference model")
    modelsComparisonMetrics.set_reaction_sets(reaction_sets=reaction_sets)

    modelsComparisonMetrics.draw_dot_plot(jaccard_distances=jaccard_distances_reactions,
                                          ratios=ratios_reactions,
                                          output_path=reactions_dotplot_path)


if __name__ == "__main__":

    print("Generating results for Bordetella pertussis")
    resultsGeneration(organism="Bpertussis",
                      genes_dotplot_path="results/Bpertussis-genes-dot-plot.jpeg",
                      reactions_dotplot_path="results/Bpertussis-reactions-dot-plot.jpeg")

    # print("Generating results for Chlorella vulgaris")
    # resultsGeneration(organism="Cvulgaris",
    #                   genes_dotplot_path="results/Cvulgaris-genes-dot-plot.jpeg",
    #                   reactions_dotplot_path="results/Cvulgaris-reactions-dot-plot.jpeg")

    print("Generating results for Lactobacillus plantarum")
    resultsGeneration(organism="Lplantarum",
                      genes_dotplot_path="results/Lplantarum-genes-dot-plot.jpeg",
                      reactions_dotplot_path="results/Lplantarum-reactions-dot-plot.jpeg")

    # print("Generating results for Toxoplasma gondii")
    # resultsGeneration(organism="Tgondii ",
    #                   genes_dotplot_path="results/Tgondii-genes-dot-plot.jpeg",
    #                   reactions_dotplot_path="results/Tgondii-reactions-dot-plot.jpeg")
