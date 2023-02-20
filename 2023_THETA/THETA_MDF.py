#%%
import numpy as np
import pandas as pd
import equilibrator_api
from equilibrator_api import ComponentContribution, Q_, ureg, Reaction
import equilibrator_pathway
from equilibrator_pathway import ThermodynamicModel
import matplotlib.pyplot as plt

print('equlibrator_api version:', equilibrator_api.__version__)
print('equlibrator_pathway version:', equilibrator_pathway.__version__)

ureg.default_format = ".2f~P"
plt.rc('axes', axisbelow=True)
ureg.setup_matplotlib(True)

comp_contrib = ComponentContribution()

# %%
comp_contrib.p_h = Q_(7)
comp_contrib.ionic_strength = Q_('250 mM')
comp_contrib.p_mg = Q_(3)

# %%
THETA = ThermodynamicModel.from_sbtab("THETA.tsv", comp_contrib=comp_contrib) 
print(THETA.net_reaction_formula)

# %%
THETA.update_standard_dgs()
mdf_result_THETA = THETA.mdf_analysis()

mdf_result_THETA.reaction_df

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
mdf_result_THETA.plot_driving_forces(ax)
ax.grid('on')
fig.savefig('mdf_result_THETA_ambient.eps')
#%%
fig, ax = plt.subplots(1, 1, figsize=(10, 15))
mdf_result_THETA.plot_concentrations(ax)
