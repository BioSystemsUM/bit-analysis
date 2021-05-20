import os
from typing import List

import pandas as pd
import subprocess
import glob
import numpy as np

from progressbar import ProgressBar
from pca_analysis import ModelAnalysis, read_models


def run_recognizer():
    for file in glob.glob('bit-analysis/comparative_func_analysis/faas/*.faa'):
        output_folder = file.split('/')[-1].split('.faa')[0]
        command = 'recognizer.py -f {} -o bit_analysis/recognizer_outs/{} -rd resources_directory --remove-spaces -dbs COG KOG'.format(
            file, output_folder)
        print(command)
        subprocess.run(command.split())

    cog_protein_descriptions = dict()

    for file in glob.glob('bit_analysis/recognizer_outs/*/COG_report.tsv'):
        name = file.split('/')[2];
        print(name)
        data = pd.read_csv(f'bit_analysis/recognizer_outs/{name}/COG_report.tsv', sep='\t')
        data = data[data['evalue'] < 10e-10]
        cog_protein_descriptions[name] = data[data['COG general functional category'] == 'METABOLISM']['DB ID'].tolist()

    return cog_protein_descriptions


def boolean_matrix(dictionary):
    cpd = list()
    for v in dictionary.values():
        cpd += v
    result = pd.DataFrame(index=range(len(set(cpd))))
    for name in list(dictionary.keys()):
        result = pd.concat([result, pd.Series([0] * len(set(cpd)), name=name)], axis=1)
    result.index = set(cpd)
    result = result.transpose()
    for k, v in dictionary.items():
        print(k)
        pbar = ProgressBar()
        for prot in pbar(v):
            result.loc[k][prot] = 1
    return result


def distance_between_rows(result, row1, row2):
    return (result.loc[row1] - result.loc[row2]).abs().sum()


def distance_matrix(df):
    result = pd.DataFrame(index=df.index, columns=df.index)
    for taxon1 in df.index:
        for taxon2 in df.index:
            result.loc[taxon1, taxon2] = distance_between_rows(df, taxon1, taxon2)
    return result


def genomes_cog_analysis():
    cpd = run_recognizer()
    b_matrix = boolean_matrix(cpd)
    b_matrix.to_csv('bit_analysis/boolean_matrix.tsv', sep='\t')
    d_matrix = distance_matrix(b_matrix)
    d_matrix.to_csv('bit_analysis/distance_matrix.tsv', sep='\t')


def models_cog_analysis():
    id2species = {'Mtub': 'Mycobacterium_tuberculosis',
                  'Sthe': 'Streptococcus_thermophilus',
                  'Xfas': 'Xylella_fastidiosa'}

    id2model = {
        'Mtub': pd.read_csv('bit_analysis/recognizer_outs/Mycobacterium_tuberculosis/COG_report.tsv', sep='\t'),
        'Sthe': pd.read_csv('bit_analysis/recognizer_outs/Streptococcus_thermophilus/COG_report.tsv', sep='\t'),
        'Xfas': pd.read_csv('bit_analysis/recognizer_outs/Xylella_fastidiosa/COG_report.tsv', sep='\t')
    }

    refseq2cog = dict()

    for k, v in id2model.items():
        df = v[['qseqid', 'DB ID']]
        df.index = ['_'.join(ide.split('_')[:2]).split('.')[0] for ide in df['qseqid']]
        refseq2cog[k] = df

    mod2reac = pd.read_csv('bit-analysis/comparative_func_analysis/models_genes.tsv', sep='\t')
    mod2reac['model_genes'] = mod2reac['model_genes'].str.split(',')
    mod2reac['cogs'] = [np.nan] * len(mod2reac)

    for i in range(len(mod2reac)):
        prefix = mod2reac.iloc[i]['model_id'][:4]
        print(prefix)
        refseq2cog[prefix].drop_duplicates(inplace=True)
        mod2reac.loc[mod2reac.index[i], 'cogs'] = ','.join([
            refseq2cog[prefix].loc[ide]['DB ID'] for ide in mod2reac.iloc[i]['model_genes']
            if ide in refseq2cog[prefix].index])

    mod2cogs = {mod2reac.iloc[i]['model_id']: mod2reac.iloc[i]['cogs'].split(',') for i in range(len(mod2reac))}
    bmatrix = boolean_matrix(mod2cogs, mod2cogs.keys())
    bmatrix.to_csv('bit_analysis/models_bmatrix.tsv', sep='\t')


def models_genes_dataframe(models: List[ModelAnalysis]):
    index = []
    genes = []
    for model_analysis in models:
        index.append(model_analysis.model.id)

        model_genes = ','.join([gene.id for gene in model_analysis.model.genes])
        genes.append(model_genes)

    df = pd.DataFrame(data=genes,
                      index=index,
                      columns=['model_genes'])

    return df


if __name__ == '__main__':
    # models_analysis = read_models(os.path.join(os.getcwd(), 'models'))
    # models_genes_df = models_genes_dataframe(models_analysis)
    # models_genes_df.to_csv(os.path.join(os.getcwd(), 'comparative_func_analysis', 'models_genes.tsv'),
    #                        sep='\t', index_label='model_id')

    genomes_cog_analysis()
    models_cog_analysis()
