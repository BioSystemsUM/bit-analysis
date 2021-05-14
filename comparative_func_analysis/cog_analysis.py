import pandas as pd
import subprocess
import glob
from progressbar import ProgressBar

for file in glob.glob('bit_analysis/*.faa'):
    output_folder = file.split('/')[-1].split('.faa')[0]
    command = 'recognizer.py -f {} -o bit_analysis/recognizer_outs/{} -rd resources_directory --remove-spaces -dbs COG KOG'.format(file, output_folder)
    print(command)
    subprocess.run(command.split())
    
report = pd.read_excel('bit_analysis/report.xlsx')
dbs = {report.iloc[i]['name'] : report.iloc[i]['domain'] for i in range(len(report))}
cog_protein_descriptions = dict()
names = list()

pbar = ProgressBar()
for file in pbar(glob.glob('bit_analysis/recognizer_outs/*/COG_report.tsv')):
    name = file.split('/')[2]; names.append(name)
    data = pd.read_csv(f'bit_analysis/recognizer_outs/{name}/COG_report.tsv', sep='\t')  #(f'bit_analysis/recognizer_outs/{name}/{dbs[name].upper()}_report.tsv')
    cog_protein_descriptions[name] = data[data['COG general functional category'] == 'METABOLISM']['DB ID'].tolist()
    #functions = set(data[f'{dbs[name].upper()} protein description'])

cpd = list()
for v in cog_protein_descriptions.values():
    cpd += v

names = set(names)
result = pd.DataFrame(index = range(len(set(cpd))))
for name in names:
    result = pd.concat([result, pd.Series([0] * 1598, name=name)],axis=1)
result.index = set(cpd)
result = result.transpose()

for k, v in cog_protein_descriptions.items():
    print(k)
    pbar = ProgressBar()
    for prot in pbar(v):
        result.loc[k][prot] = 1

def distance_between_rows(row1, row2):
    return (result.loc[row1] - result.loc[row2]).abs().sum()

def distance_matrix(df):
    result = pd.DataFrame(index=df.index, columns=df.index)
    for taxon1 in df.index:
        for taxon2 in df.index:
            result.loc[taxon1, taxon2] = distance_between_rows(taxon1, taxon2)
    return result



