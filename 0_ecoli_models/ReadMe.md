## *E. coli* genome-scale metabolic models:

* The core model ([Orth et al., *EcoSal Plus* 2010](https://doi.org/10.1128/ecosalplus.10.2.1)).
* The most updated model *i*ML1515 ([Monk et al., *Nat Biotechnol* 2017](https://doi.org/10.1038/nbt.3956)).


## Errors/inaccuracy in the model *i*ML1515:
* Transhydrogenase (THD2pp) translocates one proton instead of two ([Bizouarn et al., *Biochim Biophys Acta* 2005](https://doi.org/10.1016/j.bbabio.2005.04.004)). This error also exists in the core model.  
* Homoserine dehydrogenase (HSDy) produces homoserine from aspartate-semialdehyde irreversibly ([He et al., *Metabolic Engineering* 2020](https://doi.org/10.1016/j.ymben.2020.03.002)).  
* Glucose dehydrogenase (GLCDpp) needs pyrroloquinoline quinone (PQQ) to work ([Adamowicz et al., *Appl Environ Microbiol* 1991](https://www.ncbi.nlm.nih.gov/pubmed/1654044)).
* Pyruvate formate lyase (PFL) is active under anaeroboic condition. The same to 2-Oxobutanoate formate lyase (OBTFL), fumarate reductase (FRD2 & FRD3).  
* The same to GarK (GLYCK2), GlxK encoded glycerate kinase (GLYCK) produces 2-phosphoglycerate (2pg), not 3-phosphoglycerate ([Zelcbuch et al., *PloS One* 2015](https://doi.org/10.1371/journal.pone.0122957)).  
* Tryptophanase (TRPAS2, tnaA) should be irreversible. Although the reaction is reversible under high pyruvate and ammonia condition ([Watanabe & Snell, *PNAS* 1972](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC426635/)), such condition is likely not existing *in vivo*.   
* Isocitrate lyase (ICL) should be reversible ([MacKintosh & Nimmo, *Biochem J* 1988](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1148809/)).  
* Gene-Protein-Reaction associations of TYRL (tyrosine lyase) and THZPSN3 (thiazole phosphate synthesis) should be exchanged to each other ([Jurgenson et al., *Annu Rev Biochem* 2009](https://doi.org/10.1146/annurev.biochem.78.072407.102340)).  
* TRPS1 = TRPS3 + TRPS2. TrpA & TrpB protein, subunits of the tryptophan synthase. Because the two reaction steps are channeled between the two subunits ([Lane & Kirschner, *Biochemistry* 1991](https://pubmed.ncbi.nlm.nih.gov/1899028/)), TRPS3 & TRPS2 may be removed from the model.  
* Ethanolamine ammonia-lyase (ETHAAL, *eutBC*) requires coenzyme B12 to work.  
* Uracil degradation pathway GPR for rutC & rutD: ([Parales & Ingraham 2010](https://journals.asm.org/doi/10.1128/JB.00573-10)) vs ([Buckner et al. 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8113886/))