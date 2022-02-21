import cobra.io
from Bio import SeqIO, Entrez
from Bio.Entrez import *
from cobra.io import read_sbml_model
import os
os.chdir("../")
Entrez.email = 'id9411@alunos.uminho.pt'

def convert_gene_ids():
    Mtuberculosis_carveme_model = read_sbml_model("../models/carveme/model_Mtuberculosis_carveme_carveme.xml")
    Mtuberculosis_all_permissive_model = read_sbml_model("../models/BIT/Permissive/model_Mtuberculosis_all_permissive.xml")
    Mtuberculosis_selected_permissive_model = read_sbml_model("../models/BIT/Permissive/model_Mtuberculosis_selected_permissive.xml")
    Mtuberculosis_random1_permissive_model = read_sbml_model("../models/BIT/Permissive/model_Mtuberculosis_random1_permissive.xml")
    path_map = {Mtuberculosis_carveme_model: "../models/carveme/model_Mtuberculosis_carveme_carveme.xml", Mtuberculosis_all_permissive_model: "../models/BIT/Permissive/model_Mtuberculosis_all_permissive.xml",
                Mtuberculosis_selected_permissive_model:"../models/BIT/Permissive/model_Mtuberculosis_selected_permissive.xml", Mtuberculosis_random1_permissive_model: "../models/BIT/Permissive/model_Mtuberculosis_random1_permissive.xml"}
    models = [Mtuberculosis_carveme_model, Mtuberculosis_all_permissive_model, Mtuberculosis_selected_permissive_model, Mtuberculosis_random1_permissive_model]
    gene_id_map = {}
    # handle = efetch(db = 'nucleotide', id = "NC_000962", retmode = 'text', rettype = 'gbwithparts')
    # records = SeqIO.read(handle, format = 'gb')
    # SeqIO.write(records, 'mtuberculosis_genome.gb', format='gb')
    records = SeqIO.read('mtuberculosis_genome.gb', format='gb')
    genes_gb = [feature for feature in records.features if feature.type == 'CDS']
    for model in models:
        for gene_gb in genes_gb:
            gene_id_map[gene_gb.qualifiers['protein_id'][0].split(".")[0]] =  gene_gb.qualifiers['locus_tag'][0]
        for gene in model.genes:
            if "Rv1175c" in gene.id:
                print()
            gene_id = gene.id
            if gene.id.count("_") > 1:
                gene_id = '_'.join(gene.id.rsplit("_")[:-1])
            if gene_id in gene_id_map.keys():
                gene.id = gene_id_map[gene_id]
            else:
                print(gene_id)
            gene.name = "G_" + gene.id
            for reaction in gene.reactions:
                if gene_id in gene_id_map.keys():
                    reaction.gene_reaction_rule = reaction.gene_reaction_rule.replace(gene_id, gene_id_map[gene_id])
                else:
                    print()
            model.repair()
        cobra.io.write_sbml_model(model, path_map[model].replace(".xml","") + "_converted.xml")

convert_gene_ids()