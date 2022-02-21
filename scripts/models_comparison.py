import pandas

from ModelsComparisonMetrics import ModelsComparisonMetrics
from utils import Utils, Type
from cobra.io.sbml import read_sbml_model


def getGeneSets(models) -> dict:
    gene_sets = {}

    for model_name in models.keys():
        gene_set = []
        model = models[model_name]
        for gene in model.genes:
            gene_set.append(gene.id)
        gene_sets[model_name] = set(gene_set)
    return gene_sets


def getReactionSets(models) -> dict:
    reaction_sets = {}
    for model_name in models.keys():
        reaction_set = []
        model = models[model_name]
        for reaction in model.reactions:
            reaction_set.append(reaction.id)
        reaction_sets[model_name] = set(reaction_set)
    return reaction_sets


def getMtuberculosis_models() -> dict:
    Mtuberculosis_carveme_model = read_sbml_model("../models/carveme/model_Mtuberculosis_carveme_carveme.xml")
    Mtuberculosis_curated_model = read_sbml_model("../models/curated_models/Mtuberculosis.xml")
    Mtuberculosis_all_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_all_restrictive.xml")
    Mtuberculosis_selected_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_selected_restrictive.xml")
    Mtuberculosis_random1_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_random1_restrictive.xml")
    Mtuberculosis_random2_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_random2_restrictive.xml")
    Mtuberculosis_random3_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_random3_restrictive.xml")
    Mtuberculosis_random4_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_random4_restrictive.xml")
    Mtuberculosis_random5_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Mtuberculosis_random5_restrictive.xml")

    Mtuberculosis_models = {"carveme model": Mtuberculosis_carveme_model,
                            "reference model": Mtuberculosis_curated_model,
                            "all model": Mtuberculosis_all_restrictive_model,
                            "selected model": Mtuberculosis_selected_restrictive_model,
                            "random1 model": Mtuberculosis_random1_restrictive_model,
                            "random2 model": Mtuberculosis_random2_restrictive_model,
                            "random3 model": Mtuberculosis_random3_restrictive_model,
                            "random4 model": Mtuberculosis_random4_restrictive_model,
                            "random5 model": Mtuberculosis_random5_restrictive_model,
                            }

    for model_name in Mtuberculosis_models:
        Mtuberculosis_models[model_name].id = model_name

    for key, model in Mtuberculosis_models.items():
        if key in ("all model", "selected model", "random1 model", "random2 model", "random3 model", "random4 model",
                   "random5 model"):
            for reaction in model.reactions:
                try:
                    reaction.id = reaction.id.split("__")[0]
                except:
                    print(key)
    return Mtuberculosis_models


def getSthermophilus_models() -> dict:
    Sthermophilus_carveme_model = read_sbml_model("../models/carveme/model_Sthermophilus_carveme_carveme.xml")
    Sthermophilus_curated_model = read_sbml_model("../models/curated_models/Sthermophilus.xml")
    Sthermophilus_all_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_all_restrictive.xml")
    Sthermophilus_selected_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_selected_restrictive.xml")
    Sthermophilus_random1_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_random1_restrictive.xml")
    Mtuberculosis_random2_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_random2_restrictive.xml")
    Mtuberculosis_random3_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_random3_restrictive.xml")
    Mtuberculosis_random4_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_random4_restrictive.xml")
    Mtuberculosis_random5_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Sthermophilus_random5_restrictive.xml")

    Sthermophilus_models = {"carveme model": Sthermophilus_carveme_model,
                            "reference model": Sthermophilus_curated_model,
                            "all model": Sthermophilus_all_restrictive_model,
                            "selected model": Sthermophilus_selected_restrictive_model,
                            "random1 model": Sthermophilus_random1_restrictive_model,
                            "random2 model": Mtuberculosis_random2_restrictive_model,
                            "random3 model": Mtuberculosis_random3_restrictive_model,
                            "random4 model": Mtuberculosis_random4_restrictive_model,
                            "random5 model": Mtuberculosis_random5_restrictive_model,
                            }

    for model_name in Sthermophilus_models:
        Sthermophilus_models[model_name].id = model_name

    for key, model in Sthermophilus_models.items():
        if key in ("all model", "selected model", "random1 model", "random2 model", "random3 model", "random4 model",
                   "random5 model"):
            for reaction in model.reactions:
                try:
                    reaction.id = reaction.id.split("__")[0]
                except:
                    print(key)
    return Sthermophilus_models


def getXfastidiosa_models() -> dict:
    Xfastidiosa_carveme_model = read_sbml_model("../models/carveme/model_Xfastidiosa_carveme_carveme.xml")
    Xfastidiosa_curated_model = read_sbml_model("../models/curated_models/Xfastidiosa.xml")
    Xfastidiosa_all_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_all_restrictive.xml")
    Xfastidiosa_selected_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_selected_restrictive.xml")
    Xfastidiosa_random1_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_random1_restrictive.xml")
    Xfastidiosa_random2_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_random2_restrictive.xml")
    Xfastidiosa_random3_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_random3_restrictive.xml")
    Xfastidiosa_random4_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_random4_restrictive.xml")
    Xfastidiosa_random5_restrictive_model = read_sbml_model(
        "../models/BIT/Restrictive/model_Xfastidiosa_random5_restrictive.xml")

    Xfastidiosa_models = {"carveme model": Xfastidiosa_carveme_model,
                          "reference model": Xfastidiosa_curated_model,
                          "all model": Xfastidiosa_all_restrictive_model,
                          "selected model": Xfastidiosa_selected_restrictive_model,
                          "random1 model": Xfastidiosa_random1_restrictive_model,
                          "random2 model": Xfastidiosa_random2_restrictive_model,
                          "random3 model": Xfastidiosa_random3_restrictive_model,
                          "random4 model": Xfastidiosa_random4_restrictive_model,
                          "random5 model": Xfastidiosa_random5_restrictive_model,
                          }

    for model_name in Xfastidiosa_models:
        Xfastidiosa_models[model_name].id = model_name

    for key, model in Xfastidiosa_models.items():
        if key in ("all model", "selected model", "random1 model", "random2 model", "random3 model", "random4 model",
                   "random5 model"):
            for reaction in model.reactions:
                try:
                    reaction.id = reaction.id.split("__")[0]
                except:
                    print(key)
    return Xfastidiosa_models


def resultsGeneration_Mtuberculosis() -> None:
    Mtuberculosis_models = getMtuberculosis_models()
    modelsComparisonMetrics = ModelsComparisonMetrics(Mtuberculosis_models)
    # modelsComparisonMetrics.draw_venn_diagram(Type.GENES, "./results/Mtuberculosis-genes-venn-diagram.jpeg")
    # modelsComparisonMetrics.draw_venn_diagram(Type.METABOLITES, "./results/Mtuberculosis-metabolites-venn-diagram.jpeg")
    # modelsComparisonMetrics.draw_venn_diagram(Type.REACTIONS, "./results/Mtuberculosis-reactions-venn-diagram.jpeg")

    # modelsComparisonMetrics.write_reaction_genes_metabolites_into_csv(organism="Mycobacterium tuberculosis",
    #                                                                   output_path_non_curated_models="./results/Mtuberculosis-stats.csv",
    #                                                                   output_path_curated_models="")

    # gene_sets = getGeneSets(models=Mtuberculosis_models)
    # modelsComparisonMetrics.set_gene_sets(gene_sets=gene_sets)

    reaction_sets = getReactionSets(models=Mtuberculosis_models)
    modelsComparisonMetrics.set_reaction_sets(reaction_sets=reaction_sets)

    jaccard_distances_genes = modelsComparisonMetrics.calculate_jaccard_distances(
        reference_model_name="reference model",
        type=Type.REACTIONS)
    ratios_genes = modelsComparisonMetrics.calculate_ratios(
        reference_model_name="reference model",
        type=Type.REACTIONS)

    # gene_sets.pop("reference model")
    # modelsComparisonMetrics.set_gene_sets(gene_sets=gene_sets)
    genes_dotplot_path = "results/Mtuberculosis-reactions-dot-plot.jpeg"
    modelsComparisonMetrics.draw_dot_plot(jaccard_distances=jaccard_distances_genes,
                                          ratios=ratios_genes,
                                          output_path=genes_dotplot_path)


def resultsGeneration_Sthermophilus() -> None:
    Sthermophilus_models = getSthermophilus_models()
    modelsComparisonMetrics = ModelsComparisonMetrics(Sthermophilus_models)
    # modelsComparisonMetrics.draw_venn_diagram(Type.GENES, "./results/Sthermophilus-genes-venn-diagram.jpeg")
    # modelsComparisonMetrics.draw_venn_diagram(Type.METABOLITES, "./results/Sthermophilus-metabolites-venn-diagram.jpeg")
    # modelsComparisonMetrics.draw_venn_diagram(Type.REACTIONS, "./results/Sthermophilus-reactions-venn-diagram.jpeg")

    # modelsComparisonMetrics.write_reaction_genes_metabolites_into_csv(organism="Streptococcus thermophilus",
    #                                                                   output_path_non_curated_models="./results/Sthermophilus-stats.csv",
    #                                                                   output_path_curated_models="")

    reaction_sets = getReactionSets(models=Sthermophilus_models)
    modelsComparisonMetrics.set_reaction_sets(reaction_sets=reaction_sets)

    jaccard_distances_genes = modelsComparisonMetrics.calculate_jaccard_distances(
        reference_model_name="reference model",
        type=Type.REACTIONS)
    ratios_genes = modelsComparisonMetrics.calculate_ratios(
        reference_model_name="reference model",
        type=Type.REACTIONS)

    genes_dotplot_path = "results/Sthermophilus-reactions-dot-plot.jpeg"
    modelsComparisonMetrics.draw_dot_plot(jaccard_distances=jaccard_distances_genes,
                                          ratios=ratios_genes,
                                          output_path=genes_dotplot_path)


def resultsGeneration_Xfastidiosa() -> None:
    Xfastidiosa_models = getXfastidiosa_models()
    modelsComparisonMetrics = ModelsComparisonMetrics(Xfastidiosa_models)
    # modelsComparisonMetrics.draw_venn_diagram(Type.GENES, "./results/Xfastidiosa-genes-venn-diagram.jpeg")
    # modelsComparisonMetrics.draw_venn_diagram(Type.METABOLITES, "./results/Xfastidiosa-metabolites-venn-diagram.jpeg")
    # modelsComparisonMetrics.draw_venn_diagram(Type.REACTIONS, "./results/Xfastidiosa-reactions-venn-diagram.jpeg")

    # modelsComparisonMetrics.write_reaction_genes_metabolites_into_csv(organism="Xylella fastidiosa",
    #                                                                   output_path_non_curated_models="./results/Xfastidiosa-stats.csv",
    #                                                                   output_path_curated_models="")

    reaction_sets = getReactionSets(models=Xfastidiosa_models)
    modelsComparisonMetrics.set_reaction_sets(reaction_sets=reaction_sets)

    jaccard_distances_genes = modelsComparisonMetrics.calculate_jaccard_distances(
        reference_model_name="reference model",
        type=Type.REACTIONS)
    ratios_genes = modelsComparisonMetrics.calculate_ratios(
        reference_model_name="reference model",
        type=Type.REACTIONS)

    genes_dotplot_path = "results/Xfastidiosa-reactions-dot-plot.jpeg"
    modelsComparisonMetrics.draw_dot_plot(jaccard_distances=jaccard_distances_genes,
                                          ratios=ratios_genes,
                                          output_path=genes_dotplot_path)


# def resultsGeneration_executionTime () -> None:
#     execution_time_values = pandas.read_excel("./results/execution-time-values.xlsx")
#     ModelsComparisonMetrics.draw_barplot(execution_time_values, "./results/execution-time-values-barplot.jpeg")

if __name__ == "__main__":
    # print("Generating results related with the execution time")
    # resultsGeneration_executionTime()

    print("Generating results for Mycobacterium tuberculosis")
    resultsGeneration_Mtuberculosis()

    print("Generating results for Streptococcus thermophilus")
    resultsGeneration_Sthermophilus()

    print("Generating results for Xylella fastidiosa")
    resultsGeneration_Xfastidiosa()
