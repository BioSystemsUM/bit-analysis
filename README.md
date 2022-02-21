# BiGG Integration Tool Analysis

Multivariate analysis of Genome-Scale Metabolic Models (GSMMs) derived from [BiGG Models](http://bigg.ucsd.edu/) database.

### **Workflow**

1. [_merlin_](https://merlin-sysbio.org/) 's [BiGG Integration Tool](https://github.com/merlin4-sysbio) 
was used to derive GSMMs from [BiGG Models](http://bigg.ucsd.edu/) database. 
BiGG Integration Tool allows one to reconstruct draft GSMMs from BiGG models using BLAST similarity searches. 
A draft reconstruction can be performed using as many BiGG template-models as you would like.
In this work, we have performed several draft reconstructions for three organisms, namely _Streptococcus thermophilus_,
   _Xylella fastidiosa_, and _Mycobacterium tuberculosis_.
The draft GSMMs have been reconstructed using the following templates: 
   - whole BiGG database; 
   - a set of functionally close models. Close models have been inferred from the COG metabolic functional annotation 
     (see bellow);
   - a set of randomly selected models (this procedure was repeated 5 times).

2. In addition, [_carveme_](https://github.com/cdanielmachado/carveme) 
   was also used to generate a GSMM for each organism using the BiGG Models database.
   
3. The multivariate analysis performed for the resulting models can be consulted in this repository.

### **Analysis**
   
#### **Comparative Genome Functional Analysis** 
1. NCBI's RefSeq database was used to retrieve the proteome sequence 
   of each organism represented in the BiGG Models database.
   
2. A COG functional annotation was performed for each organism using 
   [_reCOGnizer_](https://github.com/iquasere/reCOGnizer). 
   Then, COG metabolic functional annotation was used for further analysis.

3. Metabolic COG identifiers were used to perform a Principal Component Analysis (PCA). 
   For that, a binary matrix was constructed in which each row represented an organism of BiGG database and each column
   represented a COG identifier. All matrix cells were filled with presence (1) or absence (0) 
   of the respective COG identifier. Variables having low variance were filtered out and standardization was applied, 
   using ```VarianceThreshold``` and ```StandardScaler``` methods from [sckit-learn](), respectively.
   A ```PCA``` analysis was performed in this matrix using also [sckit-learn]().
   
4. The analysis of variance for the _reCOGnizer_ metabolic functional annotation was then compared with each organism 
   domain and phylum.
   

#### **Model Analysis** 
1. Relations among models have been investigated using Venn diagrams and PCA analysis.
   
2. [COBRApy]() package has been used to read the draft models. 
   Reactions, metabolites, and genes shared across the draft models have been inferred. 
   Furthermore, COG functional annotation of all models' genes have also been used to assess 
   the draft models' similarity. 

3. [matplotlib_venn]() package has been used to visualize shared reactions, metabolites, genes and COGs.

4. [sckit-learn]() and [matplotlib]() packages have been used to run and plot multivariate PCA analysis 
   over reactions, metabolites and metabolic COG identifiers associated to the models' genes.
   

### **Repository**
- genomes_analysis: results of the functional analysis performed to the organisms available in BiGG Models.
- model_analysis: results of the multivariate analysis performed to the draft models.
- models: draft models and some statistics 
- resources: resources for the functional analysis performed to the organisms available in BiGG Models.
- cog_analysis: module with functions to run the COG functional annotation using _reCOGnizer_
- genome_resources: module with several tools to obtain RefSeq proteome sequences for all organisms represented in BiGG Models database.
- pca_analysis: module with functions to run the PCA analysis over the models' reactions, metabolites and COG annotation
- random_models: module to select random BiGG models as template for the three organisms.
- utils: utilities
- venn_diagram: module with functions to run Venn's diagrams over the models' reactions, metabolites and COG annotation


### **Requirements**
- [merlin](https://merlin-sysbio.org/) and [BiGG Integration Tool](https://github.com/merlin4-sysbio) - draft models
- [carveme](https://github.com/cdanielmachado/carveme) - draft models
- [reCOGnizer](https://github.com/iquasere/reCOGnizer) - COG functional annotation
- [Biopython](https://github.com/biopython/biopython) - proteome sequences retrieval
- [COBRApy](https://github.com/opencobra/cobrapy) - read models
- [pandas](https://github.com/pandas-dev/pandas) - data manipulation
- [sckit-learn](https://github.com/scikit-learn/scikit-learn) - PCA analysis
- [matplotlib](https://github.com/matplotlib/matplotlib) - plot
- [matplotlib_venn](https://github.com/konstantint/matplotlib-venn) - plot


### **Installation**
Clone this repository: ``git clone https://github.com/BioSystemsUM/bit-analysis.git`` 
Consult the previous requirements list to set up all mandatory packages and tools.


### **Utilization**
BIT can be used directly on merlin (https://merlin-sysbio.org/).
Alternatively, it can be used programatically with the command:

``java -jar -Xmx5G -Xss512M -XX:+HeapDumpOnOutOfMemoryError -Djavax.xml.accessExternalDTD=all -Dworkdir=workdir/ Blast.jar " + params)``

Params:
"option"

, where option can have the following values:
1 - Use the complete BiGG database
2 - Select specific models as template
3 - Use three random models as template

"template_models"
the name of the models to use as template
 
"eValue", "bitScore", "queryCoverage", 

blast parameters for e-value, bit-score, and query coverage
 
 
"includeReactionsWithoutGPRBigg", "includeReactionsWithoutGPR"

Boolean parameters to allow including reactions with incomplete or without GPR in the BiGG database, or not.


The python script used to generate the models of this work is available in the "BIT tool files" folder - modelGenerator.py
