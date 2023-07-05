This is a repository for supplementary scripts and notebooks of publications.  
_If GitHub failed to render the Jupyter notebooks, saying 'Sorry, something went wrong. Reload?'. [nbviewer](https://nbviewer.jupyter.org/) could be used by just entering the address there._  

---
### 0. [0_ecoli_models](0_ecoli_models)
The directory contains *E. coli* genome-scale metabolic model: the **core model** and the most updated model ***i*ML1515**, and some notes on the models. 

### 1. [2020_formaldehyde condensation](2020_formaldehyde%20condensation)
The directory contains Jupyter notebooks of modelling _in vivo_ formaldehyde-THF
condensation ([He *et al.*, *Metabolites* 2020](https://doi.org/10.3390/metabo10020065)). 

The formaldehyde-THF spontaneous condensation reaction was added to the 
adjusted *E. coli* genome-scale metabolic model [*i*ML1515](https://doi.org/10.1038/nbt.3956). 
Flux balance analysis and [phenotypic phase plane](https://doi.org/10.1002/bit.10047) 
calculations were conducted using [COBRApy](https://doi.org/10.1186/1752-0509-7-74).


### 2. [2020_Promiscuous aldolases](2020_Promiscuous%20aldolases)
The directory contains Jupyter notebooks of MDF and FBA modeling of the newly designed formaldehyde assimilation pathway, the homoserine cycle ([He *et al.*, *Metabolic Engineering* 2020](https://doi.org/10.1016/j.ymben.2020.03.002)).

### 3. [2020_GED](2020_GED)
The directory contains Jupyter notebook of FBA modeling of bioproduction via the newly identified GED (Gnd-Entner-Doudoroff) cycle from CO<sub>2</sub> ([Satanowski & Dronsella *et al.*, *Nature Communications* 2020](https://doi.org/10.1038/s41467-020-19564-5)).

### 4. [2020_Energy_Auxotroph](2020_Energy_Auxotroph/LPD_FBA.ipynb)
The directory contains Jupyter notebook of FBA modeling of different NADH producing routes within a “energy-auxotroph” strain, _&Delta;lpd_ strain ([Wenk & Schann *et al.*, *Biotechnology & Bioengineering* 2020](https://doi.org/10.1002/bit.27490)).

### 5. [2020_TaCo](2020_TaCo)
The directory contains Jupyter notebooks of FBA and MDF modeling of the tartronyl-CoA (TaCo) pathway, a new-to-nature carboxylation pathway ([Scheffen *et al.*, *Nature Catalysis* 2021](https://doi.org/10.1038/s41929-020-00557-y)).

### 6. [2022_STC](2022_STC)
The directory contains Jupyter notebook of FBA modeling of the serine threonine cycle (STC), an autocatalytic C1 assimilation pathway ([Wenk *et al.*, *bioRxiv* 2022](https://www.biorxiv.org/content/10.1101/2022.09.28.509898v1)).

### 7. [2022_formate_reduction](2022_formate_reduction)  
The directory contains scripts and models for MDF and FBA modeling of different formate reduction routes ([Nattermann *et al.*, *Nature Communications* 2023](https://www.nature.com/articles/s41467-023-38072-w)). The freezed version for the publication was also deposited on [Zenodo](https://doi.org/10.5281/zenodo.7752828)

### 8. [2023_THETA](2023_THETA)  
The directory contains scripts and models for MDF analysis of the THETA cycle.   

### 9. [2023_FALD_sensors](2023_FALD_sensors)   
The directory contains scripts and models for estimation of formaldehyde dependencies of formaldehyde metabolic sensors ([Schann *et al.*, *bioRxiv* 2023](https://www.biorxiv.org/content/10.1101/2023.06.29.547045v1)).