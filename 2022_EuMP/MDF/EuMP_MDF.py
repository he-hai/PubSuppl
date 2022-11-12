# %%
import equilibrator_api
from equilibrator_api import ComponentContribution, Q_, ureg
import equilibrator_pathway
from equilibrator_pathway import ThermodynamicModel
import matplotlib.pyplot as plt

print('equlibrator_api version:', equilibrator_api.__version__)
print('equlibrator_pathway version:', equilibrator_pathway.__version__)

ureg.default_format = ".2f~P"
plt.rc('axes', axisbelow=True)
ureg.setup_matplotlib(True)

comp_contrib = ComponentContribution()

# %% [markdown]
# ## Conditions
# pH 7.5, ionic strength 0.25 M, pMg 3.  
# 
# Upper bound for FA is set to 0.5mM.

# %%
EuMP = ThermodynamicModel.from_sbtab('EuMP.tsv', comp_contrib=comp_contrib)

EuMP.update_standard_dgs()
EuMP_mdf = EuMP.mdf_analysis()

EuMP_mdf.reaction_df

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
EuMP_mdf.plot_driving_forces(ax)
ax.grid('on')
fig.savefig('mdf_result_EuMP.eps')

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 15))
EuMP_mdf.plot_concentrations(ax)