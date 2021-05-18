from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import numpy as np
import os
from pca_analysis import read_models


def get_reactions(model):
    exch = []
    for exchange in model.reactions:
        if exchange.boundary:
            exch.append(exchange.id)
            for met in exchange.metabolites:
                if met.id.endswith('_b'):
                    exch.append(exchange.id)
                    break

    reactions = [r.id for r in model.reactions if r not in exch]
    return reactions


def get_metabolites(model):
    metabolites = [m.id for m in model.metabolites if not m.id.endswith('_b') and not m.id.endswith('__e')]
    return metabolites


def find_common(reactions_a, reactions_b, reactions_c):
    abc = set(reactions_a).intersection(set(reactions_b), set(reactions_c))

    ab_all = set(reactions_a).intersection(set(reactions_b))
    ab = [r for r in ab_all if r not in abc]

    ac_all = set(reactions_a).intersection(set(reactions_c))
    ac = [r for r in ac_all if r not in abc]

    a = [r for r in reactions_a if r not in abc and r not in ab and r not in ac]

    bc_all = set(reactions_b).intersection(set(reactions_c))
    bc = [r for r in bc_all if r not in abc]

    b = [r for r in reactions_b if r not in abc and r not in ab and r not in bc]

    c = [r for r in reactions_c if r not in abc and r not in ac and r not in bc]

    return len(a), len(b), len(ab), len(c), len(ac), len(bc), len(abc)


def create_venn(values, labels, filename, title):
    venn3(subsets=values, set_labels=labels, set_colors=("orangered", "dodgerblue", "green"))
    plt.title(title)
    plt.savefig(filename)
    plt.show()


def run(models_list, analysis='Reactions'):

    if len(models_list) != 3:
        raise ValueError('3 groups are needed to perform venn\'s diagram')

    entities = []
    if analysis == "Reactions":
        for modelAnalysis in models_list:
            entities.append(get_reactions(modelAnalysis.model))
    else:
        for modelAnalysis in models_list:
            entities.append(get_metabolites(modelAnalysis.model))

    venn_values = find_common(entities[0], entities[1], entities[2])

    organisms = [model.organism for model in models_list]
    labels = []
    for orgn in organisms:
        o = orgn.split(' ')
        label = '$\it{' + o[0][0] + '. ' + o[1] + '}$'
        labels.append(label)

    if 'carveme' in models_list[0].method:
        png_name = analysis + '_' + 'carveme' + '.png'
    else:
        png_name = analysis + '_' + models_list[0].method + '_' + models_list[0].template + '.png'

    create_venn(values=venn_values, labels=labels, filename=png_name, title=analysis)


def run_random(models_list, analysis='Reactions'):

    if len(models_list) != 3:
        raise ValueError('3 groups are needed to perform venn\'s diagram')

    all_Mt_rand = []
    all_St_rand = []
    all_Xt_rand = []

    for modelAnalysis in models_list:
        if modelAnalysis.organism == 'Mycobacterium tuberculosis':
            if analysis == 'Reactions':
                entities = get_reactions(modelAnalysis.model)
            else:
                entities = get_metabolites(modelAnalysis.model)
            all_Mt_rand.append(entities)

        elif modelAnalysis.organism == 'Streptococcus thermophilus':
            if analysis == 'Reactions':
                entities = get_reactions(modelAnalysis.model)
            else:
                entities = get_metabolites(modelAnalysis.model)
            all_St_rand.append(entities)

        else:
            if analysis == 'Reactions':
                entities = get_reactions(modelAnalysis.model)
            else:
                entities = get_metabolites(modelAnalysis.model)
            all_Xt_rand.append(entities)

    all_values_rand = []
    for i in range(len(all_Mt_rand)):
        for j in range(len(all_St_rand)):
            for k in range(len(all_Xt_rand)):
                venn_values_rand = find_common(all_Mt_rand[i], all_St_rand[j], all_Xt_rand[k])
                all_values_rand.append(venn_values_rand)

    all_values_rand = np.array(all_values_rand)
    mean_values = np.average(all_values_rand, axis=0)
    mean_values = np.rint(mean_values).astype(np.int32)

    organisms = [model.organism for model in models_list]
    labels = []
    for orgn in organisms:
        o = orgn.split(' ')
        label = '$\it{' + o[0][0] + '. ' + o[1] + '}$'
        labels.append(label)

    png_name = analysis + '_' + models_list[0].method + '_' + models_list[0].template + '.png'

    create_venn(values=mean_values, labels=labels, filename=png_name, title=analysis)


if __name__ == '__main__':
    # read all models
    models = read_models(os.path.join(os.getcwd(), 'models'))

    # PERMISSIVE MODELS
    models_permissive_all = [model for model in models if model.method == 'permissive' and model.template == 'all']
    models_permissive_selected = [model for model in models if model.method == 'permissive' and
                                  model.template == 'select']

    models_permissive_random = [model for model in models if model.method == 'permissive' and
                                model.template == 'random']

    # # Reactions
    run(models_list=models_permissive_all, analysis='Reactions')  # all
    run(models_list=models_permissive_selected, analysis='Reactions')  # selected
    run_random(models_list=models_permissive_random, analysis='Reactions')  # random

    # # Metabolites
    run(models_list=models_permissive_all, analysis='Metabolites')  # all
    run(models_list=models_permissive_selected, analysis='Metabolites')  # selected
    run_random(models_list=models_permissive_random, analysis='Metabolites')  # random

    # ###########################################################################

    # RESTRICTIVE MODELS
    models_restrictive_all = [model for model in models if model.method == 'restrictive' and model.template == 'all']
    models_restrictive_selected = [model for model in models if model.method == 'restrictive' and
                                   model.template == 'select']
    models_restrictive_random = [model for model in models if model.method == 'restrictive' and
                                 model.template == 'random']

    # Reactions
    run(models_list=models_restrictive_all, analysis='Reactions')  # all
    run(models_list=models_restrictive_selected, analysis='Reactions')  # selected
    run_random(models_list=models_restrictive_random, analysis='Reactions')  # random

    # # Metabolites
    run(models_list=models_restrictive_all, analysis='Metabolites')  # all
    run(models_list=models_restrictive_selected, analysis='Metabolites')  # selected
    run_random(models_list=models_restrictive_random, analysis='Metabolites')  # random

    # ###########################################################################

    # CARVEME MODELS

    carveme = [model for model in models if 'carveme' in model.method]
    run(models_list=carveme, analysis="Reactions")
    run(models_list=carveme, analysis="Metabolites")
