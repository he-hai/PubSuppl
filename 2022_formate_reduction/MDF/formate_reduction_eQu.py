#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display
import equilibrator_api
from equilibrator_api import ComponentContribution, Q_, ureg, Reaction
import equilibrator_pathway
from equilibrator_pathway import ThermodynamicModel

print('equlibrator_api version:', equilibrator_api.__version__)
print('equlibrator_pathway version:', equilibrator_pathway.__version__)

import warnings
warnings.filterwarnings('ignore')

ureg.default_format = ".2f~P"
plt.rc('axes', axisbelow=True)
ureg.setup_matplotlib(True)

comp_contrib = ComponentContribution()

# %%
comp_contrib.p_h = Q_(7)
comp_contrib.ionic_strength = Q_('250 mM')
comp_contrib.p_mg = Q_(3)

# %%
def MDF(psw):
    print(psw)
    model = ThermodynamicModel.from_sbtab(f"{psw}.tsv", comp_contrib=comp_contrib) 
    print('Net reaction: ', model.net_reaction_formula)

    model.update_standard_dgs()
    mdf_result = model.mdf_analysis()
    mdf = mdf_result.score
    print(f'MDF: {mdf: .2f} kJ/mol')
    display(mdf_result.reaction_df)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    mdf_result.plot_driving_forces(ax)
    ax.grid('on')
    # fig.savefig(f'mdf_result_{psw}.eps')
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    mdf_result.plot_concentrations(ax)
    plt.show()

# %% [markdown]
# ## Under default conditions 
psws = ['1_THF', '1_THF_NADH', '2_CoA', '2_CoA_NADPH', '3_Pi', '3_Pi_NADH']

for psw in psws: 
    MDF(psw)

#%%
def loadModel(psw):
    print(psw)
    model = ThermodynamicModel.from_sbtab(f"{psw}.tsv", comp_contrib=comp_contrib) 
    print('Net reaction: ', model.net_reaction_formula)
    return model

def MDF_ub(model, cmp, ub):
    model.set_bounds(cid=cmp, ub=Q_(ub, ureg.mM))
    model.update_standard_dgs()
    mdf_result = model.mdf_analysis()
    return mdf_result.score

def MDF_lb(model, cmp, lb):
    model.set_bounds(cid=cmp, lb=Q_(lb, ureg.mM))
    model.update_standard_dgs()
    mdf_result = model.mdf_analysis()
    return mdf_result.score

# %% [markdown]
# ## Formate upper bound between 1 and 50 mM
mdf_for = pd.DataFrame()
mdf_for['for_ub'] = range(1,51)
for psw in psws:
    model = loadModel(psw)
    for i in mdf_for.index: 
        mdf_for.loc[i,psw]=MDF_ub(model, 'for', mdf_for['for_ub'][i])
# %%    
fig, ax = plt.subplots()
mdf_for.plot(x='for_ub', ax=ax)
plt.setp(
    ax, ylabel='MDF (kJ/mol)', ylim=(0, 22),
    xlim=(0, 50), xlabel='Formate concentration upper bound (mM)'
)
plt.savefig('formate_reduction_mdf.eps')