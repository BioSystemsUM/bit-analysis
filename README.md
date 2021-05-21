# bit-analysis

**Analysis of BiGG Integration Tool**

Multivariate analysis of Genome-Scale Metabolic Models (GSMMs) derived from [BiGG Models](http://bigg.ucsd.edu/) database.

**Workflow**

1. [_merlin_](https://merlin-sysbio.org/) 's [BiGG Integration Tool](https://github.com/merlin4-sysbio) 
was used to derive GSMMs from [BiGG Models](http://bigg.ucsd.edu/) database. 
BiGG Integration Tool allows one to reconstruct draft GSMMs from BiGG models using BLAST similarity searches. 
A draft reconstruction can be performed using as many BiGG template-models as you would like.
In this work, we have performed several draft reconstructions for three organisms, namely _Streptococcus thermophilus_,
   _Xylella fastidiosa_, and _Mycobacterium tuberculosis_.
The draft GSMMs have been reconstructed using the following templates: 
   - whole BiGG database; 
   - a set of functionally close models;
   - a set of randomly selected models (this procedure was repeated 5 times).

2. In addition, [_carveme_](https://github.com/cdanielmachado/carveme) 
   was also used to generate a GSMM for each organism using the BiGG Models database.
   
3. The multivariate analysis performed for the resulting models can be consulted in this repository.

**Analysis**