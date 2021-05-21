import glob
import os
import subprocess
from typing import List

import numpy as np
import pandas as pd

from pca_analysis import ModelAnalysis
from utils import read_models


def run_recognizer(work_dir):
    for file in glob.glob(f'{work_dir}/faas/*.faa'):
        output_folder = file.split('/')[-1].split('.faa')[0]
        command = f'recognizer.py -f {file} -o {work_dir}/recognizer_outs/{output_folder} -rd resources_directory --remove-spaces -dbs COG KOG'
        print(command)
        subprocess.run(command.split())

    cog_protein_descriptions = dict()

    for file in glob.glob(f'{work_dir}/recognizer_outs/*/COG_report.tsv'):
        name = file.split('/')[2]
        print(name)
        data = pd.read_csv(f'{work_dir}/recognizer_outs/{name}/COG_report.tsv', sep='\t')
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
        for prot in v:
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


def genomes_cog_analysis(work_dir):
    cpd = run_recognizer(work_dir)
    b_matrix = boolean_matrix(cpd)
    b_matrix.to_csv(f'{work_dir}/genomes_cog_analysis.tsv', sep='\t')
    d_matrix = distance_matrix(b_matrix)
    d_matrix.to_csv(f'{work_dir}/genomes_distance_matrix.tsv', sep='\t')


def models_cog_analysis(work_dir, models_dir):
    id2model = {
        'Mtuberculosis': pd.read_csv(f'{work_dir}/recognizer_outs/Mycobacterium_tuberculosis/COG_report.tsv', sep='\t'),
        'Sthermophilus': pd.read_csv(f'{work_dir}/recognizer_outs/Streptococcus_thermophilus/COG_report.tsv', sep='\t'),
        'Xfastidiosa': pd.read_csv(f'{work_dir}/recognizer_outs/Xylella_fastidiosa/COG_report.tsv', sep='\t')
    }

    refseq2cog = dict()

    for k, v in id2model.items():
        df = v[['qseqid', 'DB ID']]
        df.index = ['_'.join(ide.split('_')[:2]).split('.')[0] for ide in df['qseqid']]
        refseq2cog[k] = df

    mod2reac = pd.read_csv(f'{models_dir}/models_genes.tsv', sep='\t')
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
    bmatrix = boolean_matrix(mod2cogs)
    bmatrix.to_csv(f'{models_dir}/models_cog_analysis.tsv', sep='\t')


def models_genes_dataframe(models: List[ModelAnalysis]):
    index = []
    genes = []
    for model_analysis in models:
        index.append(model_analysis.model.id)

        if 'carveme' in model_analysis.model_id.lower():

            model_genes = ','.join([gene.id[:-2] for gene in model_analysis.model.genes])

        else:

            model_genes = ','.join([gene.id for gene in model_analysis.model.genes])

        genes.append(model_genes)

    df = pd.DataFrame(data=genes,
                      index=index,
                      columns=['model_genes'])

    return df


def write_models_genes(models_dir: str, analysis_dir: str):
    models_analysis = read_models(models_dir)
    models_genes_df = models_genes_dataframe(models_analysis)
    models_genes_df.to_csv(os.path.join(analysis_dir, 'models_genes.tsv'), sep='\t', index_label='model_id')


if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.realpath(__file__))
    # genomes_cog_analysis(f'{base_dir}/genomes_analysis')
    write_models_genes(os.path.join(base_dir, 'models'), os.path.join(base_dir, 'model_analysis', 'pca_cogs'))
    # models_cog_analysis(f'{base_dir}/genomes_analysis', f'{base_dir}/model_analysis/pca_cogs')
