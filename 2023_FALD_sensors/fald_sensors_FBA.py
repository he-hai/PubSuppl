# %%
# -*- coding: utf-8 -*-
import sys, os
import numpy as np
import pandas as pd
import cobra
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.family"] = "Arial"
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 24
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['legend.title_fontsize'] = 15
plt.rcParams['legend.fontsize'] = 13
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['savefig.bbox'] = 'tight'

print('Python version:', sys.version)
print('numpy version:', np.__version__)
print('pandas version:', pd.__version__)
print('cobrapy version:', cobra.__version__)

# %%
def KORxn(model: cobra.Model,
          rxns2KO: list):
    """Function for knocking out reactions."""
    for ID in rxns2KO:
        reaction = model.reactions.get_by_id(ID)
        reaction.knock_out()

# %%
def pfba(model: cobra.Model):
    cobra.flux_analysis.pfba(model)
    
# %%
def flux2file(model: cobra.Model, 
              psw, product, output_dir='tmp'):
    """Function of exporting flux data."""
    n = len(model.reactions)
    modelMatrix = np.empty([n, 9], dtype = object)
    for i in range(len(model.reactions)):
        x = model.reactions[i]
        modelMatrix[i, 0] = i + 1
        modelMatrix[i, 1] = x.id
        modelMatrix[i, 2] = x.name
        modelMatrix[i, 3] = x.reaction
        modelMatrix[i, 4] = x.subsystem
        modelMatrix[i, 5] = x.lower_bound
        modelMatrix[i, 6] = x.upper_bound
        modelMatrix[i, 7] = x.flux
        modelMatrix[i, 8] = abs(x.flux)
        
    df = pd.DataFrame(data = modelMatrix, 
                      columns = ['N', 'RxnID', 'RxnName', 'Reaction', 'SubSystem', 
                      'LowerBound', 'UpperBound', 'Flux-core', 'abs(Flux)'])
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    filepath = os.path.join(output_dir, '{}_{}.xlsx'.format(psw, product))
    df.to_excel(filepath, index=False)

# %% [markdown]
# ## Model background
# 
#  * Using _E. coli_ full model *i*ML1515.
#  * Modified the transhydrogenase reaction (THD2pp) from 2 protons translocation to 1.
#  * Changed HSDy (homoserine DH) to be irrversible, towards to homS.
#  * Changed ICL (isocitrate lyase) to be reversible. 
#  * Changed TRPAS2 (Tryptophanase) to be irreversible, towards tryptophane degradation. 
#  * Base deletions: PFL, OBTFL, POR5 (pyruvate synthase, Ferredoxin), FDH4pp, FDH5pp, GLYCK (glycerate kinase, to 3pg), FRD2 and FRD3.
#  * FolD is reversible, MTHFC and MTHFD, equivalent to Fch and MtdA.
#  * FTL reaction is existed, FTHFLi (irrversible)
#  * GlyA (GHMT2r) is reversible.
#  * The model has a sarcosine oxidase reaction (SARCOX, solA), but no sarcosine exchange.  
# 
# %%
# Load full model and modification 
model = cobra.io.load_json_model(r'..\0_ecoli_models\iML1515.json')

model.reactions.THD2pp.add_metabolites({"h_p": 1, "h_c": -1})
model.reactions.HSDy.bounds = (-1000, 0)
model.reactions.ICL.bounds = (-1000, 1000)
model.reactions.TRPAS2.bounds = (0, 1000)

wt = model.copy()

KORxn_base = ['POR5', 'GLYCK', 'FDH4pp', 'FDH5pp',
              'PFL', 'OBTFL', 'GART', 'DRPA', 'PAI2T',
              'FRD2','FRD3',]

KORxn(model, KORxn_base)

model.reactions.EX_glc__D_e.bounds = (0,0)

# rxn=cobra.Reaction(id='EX_sarcs',name='Sarcosine exchange',lower_bound=0,upper_bound=0)
# model.add_reaction(rxn)
# rxn.add_metabolites({'sarcs_c':-1})
# %% [markdown]
# ## RuMP sensor 
rump = model.copy()

h6p_c = cobra.Metabolite(
    'h6p_c',
    formula='C6H11O9P',
    name='arabino-3-hexulose 6-phosphate',
    compartment='c'
)

rxn=cobra.Reaction(id='HPS',name='3-hexulose-6-phosphate synthase',lower_bound=-1000,upper_bound=1000)
rump.add_reaction(rxn)
rxn.add_metabolites({'ru5p__D_c':-1,'fald_c':-1,h6p_c:1})

rxn=cobra.Reaction(id='PHI',name='6-phospho-3-hexuloisomerase',lower_bound=-1000,upper_bound=1000)
rump.add_reaction(rxn)
rxn.add_metabolites({h6p_c:-1,'f6p_c':1})

KO_rump = [
    'FALDH2',  # frmA
    'TKT1','TKT2',
    'FBP',  # fbp, glpX
    'G6PDH2r',  # zwf
]
KORxn(rump,KO_rump)

rump.reactions.EX_xyl__D_e.lower_bound = -1000
rump.reactions.EX_succ_e.lower_bound = -1000

# E4P source 
rxn=cobra.Reaction(id='EX_E4P',name='EX E4P',lower_bound=-1000,upper_bound=0)
rump.add_reaction(rxn)
rxn.add_metabolites({'e4p_c':-1})

# %%
# check growth without sarcosine/formaldehyde first 
with rump as m: 
    # need to KO F6PA: dha + g3p <=> f6 to create the sensor strain 
    m.reactions.F6PA.knock_out()
    pfba(m)
    print(m.summary())
    # flux2file(m,'rump','0')
# %%
with rump as m: 
    # m.reactions.EX_sarcs.lower_bound = -1
    m.reactions.EX_fald_e.lower_bound = -1
    m.reactions.F6PA.knock_out()
    print(m.summary())
    # flux2file(m,'rump','1')

# %% [markdown]
# When formaldehyde uptake is the constraint to the growth, 
# the slope between the minimum required formaldehyde uptake rate (y-axis) and growth rate (x-axis) 
# is the formaldehyde requirements/dependency of the model, i.e., 
# $$\dfrac{mmol/gCDW/h}{1/h} = \dfrac{mmol}{gCDW}$$

# %%
data = pd.DataFrame(
    columns=['Growth rate','FALD uptake'],
    dtype=float,
)
data['Growth rate'] = np.arange(0.5,step=0.1)

def req_calc(
        model: cobra.Model,
        df: pd.DataFrame,
        strain: str
):
    model.reactions.EX_fald_e.lower_bound = -1000
    model.objective ={model.reactions.EX_fald_e:-1}
    model.objective_direction='min'

    for i in df.index: 
        gr = df.loc[i,'Growth rate']
        # print(gr)
        model.reactions.BIOMASS_Ec_iML1515_core_75p37M.bounds = (gr,gr)
        pfba(model)
        # print(model.summary())
        df.loc[i,'FALD uptake'] = abs(model.reactions.EX_fald_e.flux)

    slope = (df.iloc[4,1]-df.iloc[0,1])/(df.iloc[4,0]-df.iloc[0,0])

    fig, ax = plt.subplots(figsize=(5,5))
    sns.lineplot(
        x='Growth rate',y='FALD uptake',data=df,
        marker='o',ax=ax,
    )
    ax.fill_between(
        x=df['Growth rate'],y1=0.5,y2=df['FALD uptake'],
        facecolor='b',alpha=0.5,
    )
    plt.setp(
        ax,xlim=(0,0.4),ylim=(0,0.5),
        xlabel='Growth rate (1/h)',
        ylabel='FALD uptake rate (mmol/gCDW/h)',
        title=f'{strain}\n slope={slope:.3f} mmol/gCDW',
    )
    plt.savefig(f'{strain} FALD dependecy.pdf')

# %%
with rump as m:
    m.reactions.F6PA.knock_out()
    rump_df = data.copy()
    req_calc(m,rump_df,'RuMP')

# %% [markdown]
# ## LtaE sensor 
ltaE = model.copy()

rxn=cobra.Reaction(id='SAL',name='serine aldolase',lower_bound=-1000,upper_bound=1000)
ltaE.add_reaction(rxn)
rxn.add_metabolites({'gly_c':-1,'fald_c':-1,'ser__L_c':1})

KO_ltaE = [
    'FALDH2',  # frmA
    'AHGDx','PGCD', # serA
    'GHMT2r','THFAT'  # glyA
]
KORxn(ltaE,KO_ltaE)

ltaE.reactions.EX_glc__D_e.lower_bound = -1000
ltaE.reactions.EX_gly_e.lower_bound = -1000
# %%
# check growth without sarcosine/formaldehyde first 
with ltaE as m: 
    pfba(m)
    print(m.summary())
    # flux2file(m,'ltaE','0')
# %%
with ltaE as m: 
    m.reactions.EX_fald_e.lower_bound = -1
    print(m.summary())
    # flux2file(m,'ltaE','1')

# %%
with ltaE as m:
    ltaE_df = data.copy()
    req_calc(m,ltaE_df,'LtaE')

# %% [markdown]
# ## HAL sensor 
HAL = model.copy()

hob_c = cobra.Metabolite(
    'hob_c',
    formula='C4H5O4',
    name='4-hydroxy-2-oxobutanoate',
    compartment='c'
)

rxn=cobra.Reaction(id='HAL',name='HOB aldolase',lower_bound=-1000,upper_bound=1000)
HAL.add_reaction(rxn)
rxn.add_metabolites({'pyr_c':-1,'fald_c':-1,hob_c:1})

rxn=cobra.Reaction(id='HAT',name='HOB aminotransferase',lower_bound=-1000,upper_bound=1000)
HAL.add_reaction(rxn)
rxn.add_metabolites({hob_c:-1,'glu__L_c':-1,'hom__L_c':1,'akg_c':1})

KO_hal = [
    'FALDH2',  # frmA
    'ASAD', # asd
]
KORxn(HAL,KO_hal)

HAL.reactions.EX_glc__D_e.lower_bound = -1000
HAL.reactions.EX_26dap__M_e.lower_bound = -1000
# %%
# check growth without sarcosine/formaldehyde first 
with HAL as m: 
    pfba(m)
    print(m.summary())
    # flux2file(m,'HAL','0')

# %%
with HAL as m: 
    m.reactions.EX_fald_e.lower_bound = -1
    print(m.summary())
    # flux2file(m,'HAL','1')

# %%
with HAL as m: 
    HAL_df = data.copy()
    req_calc(m,HAL_df,'HAL')

# %% [markdown]
# ## HAL with Thr
with HAL as m: 
    m.reactions.EX_thr__L_e.lower_bound = -1000
    pfba(m)
    print(m.summary())
    # flux2file(m,'HAL_Thr','0')
# %%
with HAL as m: 
    m.reactions.EX_thr__L_e.lower_bound = -1000
    m.reactions.EX_fald_e.lower_bound = -1
    print(m.summary())
    # flux2file(m,'HAL_Thr','1')

# %%
with HAL as m: 
    m.reactions.EX_thr__L_e.lower_bound = -1000
    HAL2_df = data.copy()
    req_calc(m,HAL2_df,'HAL thr')

# %%
