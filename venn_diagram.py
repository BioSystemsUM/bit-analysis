from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn2
import numpy as np
import os
import pandas as pd
from utils import read_models, parse_reaction, parse_metabolite


# Get reactions from models, excluding exchanges
def get_reactions(model):
    reactions = []
    for r in model.reactions:
        parse = parse_reaction(r)
        if parse:
            reactions.append(parse)
    return reactions


# Get metabolites from models, excluding extracellular
def get_metabolites(model):
    metabolites = []
    for m in model.metabolites:
        parse = parse_metabolite(m)
        if parse:
            metabolites.append(parse)
    return metabolites


# Get genes from models
def get_genes(model):
    genes = []
    for gene in model.genes:
        if gene.id[-2] == '_':
            genes.append(gene.id[:-2])
        else:
            genes.append(gene.id)
    return genes


# Find common reactions/metabolites/genes between models
def find_common(entities_a, entities_b, entities_c):
    abc = set(entities_a).intersection(set(entities_b), set(entities_c))

    ab_all = set(entities_a).intersection(set(entities_b))
    ab = [r for r in ab_all if r not in abc]

    ac_all = set(entities_a).intersection(set(entities_c))
    ac = [r for r in ac_all if r not in abc]

    a = [r for r in entities_a if r not in abc and r not in ab and r not in ac]

    bc_all = set(entities_b).intersection(set(entities_c))
    bc = [r for r in bc_all if r not in abc]

    b = [r for r in entities_b if r not in abc and r not in ab and r not in bc]

    c = [r for r in entities_c if r not in abc and r not in ac and r not in bc]

    return len(a), len(b), len(ab), len(c), len(ac), len(bc), len(abc)


# Create a venn's diagram for 2 or 3 groups
def create_venn(values, labels, filename, title):
    if len(labels) == 3:
        venn3(subsets=values, set_labels=labels, set_colors=("orangered", "dodgerblue", "green"))
    else:
        venn2(subsets=values, set_labels=labels, set_colors=("orangered", "dodgerblue"))
    plt.title(title)
    plt.savefig(filename)
    plt.show()


# Group
def run(models_list, analysis='Reactions', group='organism'):
    entities = []
    if len(models_list) <= 3:
        if analysis == "Reactions":
            for modelAnalysis in models_list:
                entities.append(get_reactions(modelAnalysis.model))
        elif analysis == "Metabolites":
            for modelAnalysis in models_list:
                entities.append(get_metabolites(modelAnalysis.model))
        else:
            for modelAnalysis in models_list:
                entities.append(get_genes(modelAnalysis.model))

        if len(models_list) == 3:
            venn_values = find_common(entities[0], entities[1], entities[2])
        elif len(models_list) == 2:
            venn_values = find_common(entities[0], entities[1], [])
        else:
            raise ValueError('At least 2 groups are needed to perform Venn Diagrams')

    else:
        venn_values = run_random(models_list, analysis, group)

    labels = []
    if group == 'organism':
        for model in models_list:
            org = model.organism.split(' ')
            label = '$\it{' + org[0][0] + '. ' + org[1] + '}$'
            if label not in labels:
                labels.append(label)
        png_name = 'model_analysis/venn/' + analysis + '_' + models_list[0].method + '_' \
                   + models_list[0].template + '.png'

    elif group == 'template':
        for model in models_list:
            if model.template not in labels:
                labels.append(model.template)
        png_name = 'model_analysis/venn/' + models_list[0].organism_id + '_' + analysis + '_' \
                   + models_list[0].method + '.png'

    else:
        for model in models_list:
            if model.method == 'carveme':
                labels.append(model.method)
            else:
                labels.append('bit_' + model.template)
        png_name = 'model_analysis/venn/' + 'bitVScarveme' + '_' + models_list[0].organism_id + '_' + analysis + '_' + \
                   models_list[1].method + '.png'

    create_venn(values=venn_values, labels=labels, filename=png_name, title=analysis)


# Groups models when models_list have more than 3 models (for random ones)
def run_random(models_list, analysis='Reactions', group='organism'):
    groups = defaultdict(list)
    if group == 'organism':
        for model in models_list:
            groups[model.organism_id].append(model)
    elif group == 'template':
        for model in models_list:
            groups[model.template].append(model)
    groups_list = list(groups.values())

    all_entities = []
    for g in groups_list:
        if analysis == 'Reactions':
            entities = [get_reactions(modelAnalysis.model) for modelAnalysis in g]
        elif analysis == 'Metabolites':
            entities = [get_metabolites(modelAnalysis.model) for modelAnalysis in g]
        else:
            entities = [get_genes(modelAnalysis.model) for modelAnalysis in g]
        all_entities.append(entities)

    all_values = []
    for i in range(len(all_entities[0])):
        for j in range(len(all_entities[1])):
            for k in range(len(all_entities[2])):
                venn_values = find_common(all_entities[0][i], all_entities[1][j], all_entities[2][k])
                all_values.append(venn_values)

    all_values_np = np.array(all_values)
    mean_values = np.average(all_values_np, axis=0)
    mean_values = np.rint(mean_values).astype(np.int32)

    return mean_values


# Create Venn's Diagram for COGs analysis
def run_cogs_genomes(dataset):
    df = pd.read_csv(dataset, sep='\t', header=0)
    organisms = df['organism'].tolist()
    cogs = df['cogs'].tolist()
    cogs = [cog.split(',') for cog in cogs]
    common = find_common(cogs[0], cogs[1], cogs[2])
    labels = []
    for organism in organisms:
        org = organism.split('_')
        label = '$\it{' + org[0][0] + '. ' + org[1] + '}$'
        labels.append(label)
    create_venn(values=common, labels=labels, filename='model_analysis/venn/COGs_org.png', title='COGs')


# Create Venn's Diagrams for all permissive models
def run_permissive(models_list, analysis='Reactions'):
    models_permissive_all = [model for model in models_list if model.method == 'permissive' and model.template == 'all']
    models_permissive_selected = [model for model in models_list if model.method == 'permissive' and
                                  model.template == 'select']
    models_permissive_random = [model for model in models_list if model.method == 'permissive' and
                                model.template == 'random']

    run(models_list=models_permissive_all, analysis=analysis, group='organism')  # all
    run(models_list=models_permissive_selected, analysis=analysis, group='organism')  # selected
    run(models_list=models_permissive_random, analysis=analysis, group='organism')  # random


# Create Venn's Diagrams for all restrictive models
def run_restrictive(models_list, analysis='Reactions'):
    models_restrictive_all = [model for model in models_list if
                              model.method == 'restrictive' and model.template == 'all']
    models_restrictive_selected = [model for model in models_list if model.method == 'restrictive' and
                                   model.template == 'select']
    models_restrictive_random = [model for model in models_list if model.method == 'restrictive' and
                                 model.template == 'random']

    run(models_list=models_restrictive_all, analysis=analysis, group='organism')  # all
    run(models_list=models_restrictive_selected, analysis=analysis, group='organism')  # selected
    run(models_list=models_restrictive_random, analysis=analysis, group='organism')  # random


# Create Venn's Diagrams for all models of an organism
def run_organism(models_list, organism, method):
    group = [model for model in models_list if organism in model.organism_id
             and model.method == method]
    run(models_list=group, analysis='Reactions', group='template')
    # run(models_list=group, analysis='Metabolites', group='template')
    # run(models_list=group, analysis='Genes', group='template')


# Compares CarveMe and bit selected models
def compare_carveme_bit(models_list, organism, method):
    group = []
    for model in models_list:
        if organism in model.organism_id and model.method == 'carveme':
            group.append(model)
        if organism in model.organism_id and model.method == method and model.template == 'select':
            group.append(model)

    run(models_list=group, analysis='Reactions', group='method')
    run(models_list=group, analysis='Metabolites', group='method')
    run(models_list=group, analysis='Genes', group='method')


def models_statistics(models_list):
    columns = ['method', 'template', 'organism', 'genes', 'reactions', 'exchanges', 'metabolites',
               'extracellularMetabolites']

    data = []

    for modelAnalysis in models_list:
        method = modelAnalysis.method
        organism = modelAnalysis.organism_id
        template = modelAnalysis.template
        genes = len(get_genes(modelAnalysis.model))
        reactions = len(get_reactions(modelAnalysis.model))
        exchanges_n = len(modelAnalysis.model.exchanges)
        metabolites = len(get_metabolites(modelAnalysis.model))
        met_extra_n = len([m for m in modelAnalysis.model.metabolites if m.id.endswith('_b')])

        line = [method, template, organism, genes, reactions, exchanges_n, metabolites, met_extra_n]
        data.append(line)

    df = pd.DataFrame(data, columns=columns)
    df.to_csv('model_analysis/models_statistics.tsv', sep='\t')


def run_cogs_models(dataset, method='restrictive', template='selected',
                    filename='model_analysis/venn/COGs_select_models.png'):
    df = pd.read_csv(dataset, sep='\t', header=0)
    subdf = df.loc[df['model_id'].str.contains(template) & df['model_id'].str.contains(method)]
    subdf = subdf.where(subdf != 1, subdf.columns.to_series(), axis=1)
    model_ids = subdf['model_id'].to_list()
    cogs = []
    for ind, row in subdf.iterrows():
        cogs.append(row.to_list()[1:])
    cogs = [list(filter(lambda a: a != 0, cog)) for cog in cogs]

    common = find_common(cogs[0], cogs[1], cogs[2])

    labels = []
    for name in model_ids:
        org = name.split('_')[0]
        label = '$\it{' + org[0] + '. ' + org[1:] + '}$'
        labels.append(label)

    create_venn(values=common, labels=labels, filename=filename, title='COGs')


if __name__ == '__main__':
    # # read all models
    # models = read_models(os.path.join(os.getcwd(), 'models'))
    #
    # # COMPARISON BETWEEN ORGANISMS
    #
    # # PERMISSIVE MODELS
    # run_permissive(models_list=models, analysis='Reactions')
    # run_permissive(models_list=models, analysis='Metabolites')
    #
    # # RESTRICTIVE MODELS
    # run_restrictive(models_list=models, analysis='Reactions')
    # run_restrictive(models_list=models, analysis='Metabolites')
    #
    # # CARVEME MODELS
    # carveme = [model for model in models if 'carveme' in model.method]
    #
    # run(models_list=carveme, analysis="Reactions", group='organism')
    # run(models_list=carveme, analysis="Metabolites", group='organism')
    #
    # # COMPARISON FOR EACH ORGANISM
    #
    # # M. tuberculosis
    # run_organism(models_list=models, organism='Mtub', method='permissive')
    # run_organism(models_list=models, organism='Mtub', method='restrictive')
    #
    # # S. thermophilus
    # run_organism(models_list=models, organism='Sthe', method='permissive')
    # run_organism(models_list=models, organism='Sthe', method='restrictive')
    #
    # # X. fastidiosa
    # run_organism(models_list=models, organism='Xfas', method='permissive')
    # run_organism(models_list=models, organism='Xfas', method='restrictive')
    #
    # # RUN FOR COGS
    # run_cogs_genomes('genomes_analysis/protagonists2cogs.tsv')
    run_cogs_models('model_analysis/pca_cogs/models_cog_analysis.tsv',
                    method='restrictive',
                    template='selected',
                    filename='model_analysis/venn/COGs_select_models.png')
    run_cogs_models('model_analysis/pca_cogs/models_cog_analysis.tsv',
                    method='carveme',
                    template='carveme',
                    filename='model_analysis/venn/COGs_carveme_models.png')
    #
    # # COMPARISON bit vs CARVE ME
    # compare_carveme_bit(models_list=models, organism='Mtub', method='permissive')
    # compare_carveme_bit(models_list=models, organism='Sthe', method='permissive')
    # compare_carveme_bit(models_list=models, organism='Xfas', method='permissive')
    #
    # compare_carveme_bit(models_list=models, organism='Mtub', method='restrictive')
    # compare_carveme_bit(models_list=models, organism='Sthe', method='restrictive')
    # compare_carveme_bit(models_list=models, organism='Xfas', method='restrictive')
    #
    # # create table
    # models_statistics(models_list=models)
