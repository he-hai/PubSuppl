# %%
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import cobra
plt.rc('axes', axisbelow=True)

print('Python version:', sys.version)
print('numpy version:', np.__version__)
print('pandas version:', pd.__version__)
print('cobrapy version:', cobra.__version__)

# %%
def AddRxn(model, newRxnFile):
    """Function of adding new reactions to the model."""
    n1 = len(model.reactions)
    AllAddRxn = pd.read_csv(newRxnFile, sep=',', index_col='RxnID', skipinitialspace=True)
    n2 = len(AllAddRxn)
    for i in range(n2):
        ID = AllAddRxn.index.values[i]
        addRxn = cobra.Reaction(ID)
        model.add_reactions([addRxn])
        addRxnInf = model.reactions[n1 + i]
        addRxnInf.name = AllAddRxn.loc[ID, 'RxnName']
        addRxnInf.reaction = AllAddRxn.loc[ID, 'RxnFormula']
        addRxnInf.subsystem = AllAddRxn.loc[ID, 'Subsystem']
        addRxnInf.lower_bound = AllAddRxn.loc[ID, 'LowerBound']
        addRxnInf.upper_bound = AllAddRxn.loc[ID, 'UpperBound']
    return model


# %%
def flux2file(model, product, psw, output_dir='tmp'):
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
    filepath = os.path.join(output_dir, '{}_{}.xlsx'.format(product, psw))
    df.to_excel(filepath, index=False)
    

# %% [markdown]
# Construct a cobora model  
# The model has a constraint of producing 1 of 3pg (sink_3pg), 
# set the objective to minimize ATP cost. 
model = cobra.Model()

AddRxn(model,'CBBrxns.csv')

model.objective = {model.reactions.DM_atp: 1}
model.objective_direction = 'min'
# %% [markdown]
# Set the carboxylation and oxygenation reaction ratio 
# of RuBisCO to be 3:1
rubisco_flux = model.problem.Constraint(
    model.reactions.RBPC.flux_expression - 3 * model.reactions.RBPO.flux_expression,
    lb = 0, 
    ub = 0
)
model.add_cons_vars(rubisco_flux)

# %% [markdown]
# Routes to compare (see [Scheffen et al 2021 Nat Catal](https://doi.org/10.1038/s41929-020-00557-y),
# [Trudeau et al 2018 PNAS](https://doi.org/10.1073/pnas.1812605115) & this study): 
# - Natural photorespiration (NPR)
# - Glycerate route (GLC)
# - GLycolate oxidation pathway (OX)
# - Arabinose 5-phosphate shunt (A5P) 
# - Tartronyl-CoA pathway (TACO)
# - CORE with an acetoacetyl-CoA synthetase (CORE_Lig)
# - CORE with an formate-acetoacetate CoA-transferase (CORE_CoAT)
photores = {
    'NPR': 'a_NPRrxns.csv',    # natural photorespiration
    'GLC': 'b_GLCrxns.csv',    # glycerate bypass
    'OX': 'c_OXrxns.csv',      # glycolate oxidation pathway
    'A5P': 'd_A5Prxns.csv',    # arabinose-5-phosphate shunt
    '3OHP': 'e_3OHPrxns.csv',  # 3-hydroxypropionate bypass
    'TACO': 'f_TACOrxns.csv',  # tartronyl-CoA pathway
    'CORE_Lig': 'x1_CORE_Lig_rxns.csv',
    'CORE_CoAT': 'x2_CORE_CoAT_rxns.csv',
}

cost_df = pd.DataFrame()

# %%
for psw, rxns in photores.items():
    with model as m:
        AddRxn(m, rxns)
        m.optimize()
        flux2file(m,'3pg',psw,'output')
        for cost in ['DM_atp', 'DM_e', 'Fdr', 'EX_co2', 'RBPC', 'RBPO', 'sink_3pg']:
            cost_df.loc[cost, psw] = abs(m.reactions.get_by_id(cost).flux)
    

# %% [markdown]
# Use the GCC M5 variant that hydrolyse 3.9 ATP per carboxylation 
# (Fig. 2c of [TACO](https://doi.org/10.1038/s41929-020-00557-y))

with model as m:
    AddRxn(m, 'f_TACOrxns.csv')
    m.reactions.GCC.add_metabolites({'atp': -2.9, 'adp': 2.9, 'pi': 2.9})
    print(m.reactions.GCC)
    m.optimize()
    flux2file(m,'3pg', 'TACO_M5','output')
    for cost in ['DM_atp', 'DM_e', 'Fdr', 'EX_co2', 'RBPC', 'RBPO', 'sink_3pg']:
        cost_df.loc[cost, 'TACO_M5'] = abs(m.reactions.get_by_id(cost).flux)

# %% [markdown]
# Use the GCC M5 L100N variant that hydrolyse 1.7 ATP per carboxylation 
# ([Marchal et al 2023 ACS Synth Biol](https://doi.org/10.1021/acssynbio.3c00403))

with model as m:
    AddRxn(m, 'f_TACOrxns.csv')
    m.reactions.GCC.add_metabolites({'atp': -0.7, 'adp': 0.7, 'pi': 0.7})
    print(m.reactions.GCC)
    m.optimize()
    flux2file(m,'3pg', 'TACO_L100N','output')
    for cost in ['DM_atp', 'DM_e', 'Fdr', 'EX_co2', 'RBPC', 'RBPO', 'sink_3pg']:
        cost_df.loc[cost, 'TACO_L100N'] = abs(m.reactions.get_by_id(cost).flux)
# %%
cost_df
# %% [markdown]
# Assuming 1 NAD(P)H = 2 ferredoxin,  
# and along ETC: 1 NAD(P)H = 2.5 ATP
cost_df.loc['reducing_eq',:]=cost_df.loc['DM_atp',:]/2.5 + cost_df.loc['DM_e'] + cost_df.loc['Fdr']/2
cost_df.loc['Yield',:]=1/cost_df.loc['reducing_eq',:]  # per reducing equivalent  
cost_df.loc['Rel_yield',:]=cost_df.loc['Yield',:]/cost_df.loc['Yield','NPR']*100
cost_df.to_excel('CORE costs comparison.xlsx')
cost_df.round(3)
# %%
fig, ax = plt.subplots()
cost_df.loc['Rel_yield'].plot(kind='bar', ax=ax)
ax.grid(axis='y')
plt.savefig('CORE cost comparison.eps',dpi=300,format='pdf')
plt.show()
# %%
