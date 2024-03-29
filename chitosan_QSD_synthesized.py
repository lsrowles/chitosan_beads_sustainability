#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: l.s. rowles
    
    
    
"""


import numpy as np
import pandas as pd
import lhs
import math
import scipy
from scipy import stats
from setup import setup_data
import copy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


#%% LCC/LCA modeling inputs

# outputs of interest
output_perc_mid = 50
output_perc_low = 5
output_perc_high = 95

# general parameters - import spreadsheet tabs as dataframes
general_assumptions = pd.read_excel('./assumptions_chitosan.xlsx', sheet_name = 'General', index_col = 'Parameter')
design_assumptions = pd.read_excel('./assumptions_chitosan.xlsx', sheet_name = 'Design', index_col='Parameter')
LCA_assumptions = pd.read_excel('./assumptions_chitosan.xlsx', sheet_name = 'LCA', index_col='Parameter')

# number of Monte Carlo runs
n_samples = int(general_assumptions.loc['n_samples','expected'])

expected_lifetime = int(general_assumptions.loc['expected_lifetime','expected'])
zero_LCA = np.reshape(LCA_assumptions.loc['zero_LCA',:].iloc[2:].to_numpy(dtype=float), (1,-1))

# create empty datasets to eventually store data for sensitivity analysis (Spearman's coefficients)
correlation_distributions = np.full((n_samples, n_samples), np.nan)
correlation_parameters = np.full((n_samples, 1), np.nan)
correlation_parameters = correlation_parameters.tolist()


# loading in all of the data 
result = setup_data([general_assumptions,design_assumptions,LCA_assumptions],correlation_distributions,correlation_parameters,n_samples)

# creates variables for each of the variables in the excel file

for key in result.keys():
    exec(f'{key} = result[key]')

# LCA inputs (includes factors for all TRACI impact categories)
electricity_med_voltage_impacts = np.reshape(LCA_assumptions.loc['electricity_med_voltage_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
zero_LCA = np.reshape(LCA_assumptions.loc['zero_LCA',:].iloc[2:].to_numpy(dtype=float), (1,-1))
glacial_acetic_acid_impacts = np.reshape(LCA_assumptions.loc['glacial_acetic_acid_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
sodium_hydroxide_impacts = np.reshape(LCA_assumptions.loc['sodium_hydroxide_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
DI_water_impacts = np.reshape(LCA_assumptions.loc['DI_water_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
methylpentane2_impacts = np.reshape(LCA_assumptions.loc['methylpentane2_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
azodicarbonamide_impacts = np.reshape(LCA_assumptions.loc['azodicarbonamide_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
chitosan_impacts = np.reshape(LCA_assumptions.loc['chitosan_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
glutaraldehyde_impacts = np.reshape(LCA_assumptions.loc['glutaraldehyde_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
polysorbate_impacts = np.reshape(LCA_assumptions.loc['polysorbate_impacts',:].iloc[2:].to_numpy(dtype=float), (1,-1))
anionic_resin = np.reshape(LCA_assumptions.loc['anionic_resin',:].iloc[2:].to_numpy(dtype=float), (1,-1))

 
#%% general_procedure

general_material_costs = chitosan_amount * chitosan_cost / 50 \
        + glacial_acetic_acid_amount * glacial_acetic_acid_cost / 100 \
        + DI_water_syn/1000 * DI_water_cost \
        + NaOH_amount * NaOH_cost * (0.039 / 0.5) \
        + DI_water_rinsing * DI_water_cost \
        + glutaraldehyde_amount/10 * glutaraldehyde_50_cost # USD / synthesis cycle
    
general_energy_demand = hotplate_energy / 1000 * chitosan_heat_time * 0.25 # kWh / synthesis cycle, assuming only 25% of max energy demand used
general_energy_cost = general_energy_demand * electricity_cost # USD / sythesis cycle

general_material_impacts = chitosan_amount * chitosan_impacts / 1000 \
        + glacial_acetic_acid_amount * glacial_acetic_acid_impacts * 1.05 / 1000 \
        + (DI_water_syn/1000 + DI_water_rinsing) * DI_water_impacts \
        + NaOH_amount * sodium_hydroxide_impacts * 2.13 * 0.039 \
        + glutaraldehyde_amount * glutaraldehyde_impacts * 1.06 * 0.05 # impact / synthesis cycle

general_energy_impact = general_energy_demand * electricity_med_voltage_impacts # impact / sythesis cycle
#%% PCPM
PCPM_material_costs = methylpentane2_amount * methylpentane2_cost / 100 \
        + tween20_PCPM_amount * polysorbate_20_cost / 100 \
        + general_material_costs  # USD / synthesis cycle

PCPM_material_impacts = methylpentane2_amount * methylpentane2_impacts * 0.653 / 1000 \
        + tween20_PCPM_amount * polysorbate_impacts * 1.10 / 1000 \
        + general_material_impacts                                          #impact / synthesis cycle

PCPM_energy_impact = general_energy_impact                                #impact / synthesis cycle

PCPM_cost = (PCPM_material_costs + general_energy_cost) / synthesis_yield # USD / g chitosan
PCPM_impacts = (PCPM_material_impacts + general_energy_impact) / synthesis_yield # impact / g chitosan


#%% PCPA

PCPA_material_costs = azocarboxamide_amount * azodicarboxamide_97_cost / 100 \
        + general_material_costs                                            # USD / synthesis cycle

PCPA_material_impacts = azocarboxamide_amount * azodicarbonamide_impacts / 1000 \
        + general_material_impacts                                          #impact / synthesis cycle
        
PCPA_energy_impact   = general_energy_impact                                #impact / synthesis cycle


PCPA_cost = (PCPA_material_costs + general_energy_cost) / synthesis_yield # USD / g chitosan
PCPA_impacts = (PCPA_material_impacts + general_energy_impact) / synthesis_yield # impact / g chitosan
#%% PCPT

PCPT_material_costs = tween20_PCPT / 1.10 * polysorbate_20_cost / 100  \
        + general_material_costs                                             # USD / synthesis cycle
        
PCPT_material_impacts = tween20_PCPT *  polysorbate_impacts / 1000 \
        + general_material_impacts                                           #impact / synthesis cycle

PCPT_energy_impact   = general_energy_impact                                 #impact / synthesis cycle

PCPT_cost = (PCPT_material_costs + general_energy_cost) / synthesis_yield # USD / g chitosan
PCPT_impacts = (PCPT_material_impacts + general_energy_impact) / synthesis_yield # impact / g chitosan

PCPT_impacts_norm = (PCPT_impacts / (anionic_resin/1e3))
#%% outputs

costs_labels = ('PCPM_cost', 'PCPA_cost', 'PCPT_cost')
costs = (PCPM_cost, PCPA_cost, PCPT_cost)
impacts_labels = ('PCPM_impacts', 'PCPA_impacts', 'PCPT_impacts')
impacts = (PCPM_impacts, PCPA_impacts, PCPT_impacts)


#%% outputs
writer = pd.ExcelWriter('chitosan_synthesis_outputs.xlsx', engine='xlsxwriter')

df_cost = pd.DataFrame({k:v.flatten() for k,v in zip(costs_labels, costs)})

df_impact_AP = pd.DataFrame({k:v.T[0].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_GWP = pd.DataFrame({k:v.T[1].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_EFW = pd.DataFrame({k:v.T[2].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_EP = pd.DataFrame({k:v.T[3].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_HTC = pd.DataFrame({k:v.T[4].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_HTNC = pd.DataFrame({k:v.T[5].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_ODP = pd.DataFrame({k:v.T[6].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_PMPF = pd.DataFrame({k:v.T[7].flatten() for k,v in zip(impacts_labels, impacts)})
df_impact_MIR = pd.DataFrame({k:v.T[8].flatten() for k,v in zip(impacts_labels, impacts)})

df_cost.to_excel(writer, sheet_name = 'costs')

df_impact_AP.to_excel(writer, sheet_name = 'impacts_AP')

df_impact_GWP.to_excel(writer, sheet_name = 'impacts_GWP')
df_impact_EFW.to_excel(writer, sheet_name = 'impacts_EFW')
df_impact_EP.to_excel(writer, sheet_name = 'impacts_EP')
df_impact_HTC.to_excel(writer, sheet_name = 'impacts_HTC')
df_impact_HTNC.to_excel(writer, sheet_name = 'impacts_HTNC')
df_impact_ODP.to_excel(writer, sheet_name = 'impacts_ODP')
df_impact_PMPF.to_excel(writer, sheet_name = 'impacts_PMPF')
df_impact_MIR.to_excel(writer, sheet_name = 'impacts_MIR')

#%%sensitivity analysis
all_inputs_names = ('chitosan_cost', 'polysorbate_20_cost', 'azodicarboxamide_97_cost', 'methylpentane2_cost', 
              'glacial_acetic_acid_cost', 'glutaraldehyde_50_cost', 'hotplate_cost', 'hotplate_energy', 
              'NaOH_cost', 'DI_water_cost', 'chitosan_amount', 'glacial_acetic_acid_amount', 'process_temperature',
              'NaOH_amount', 'overnight_settling_time', 'DI_water_syn', 'glutaraldehyde_amount', 'DI_water_rinsing', 
              'methylpentane2_amount', 'tween20_PCPM_amount', 'azocarboxamide_amount', 'tween20_PCPT', 'chitosan_heat_time', 
              'electricity_cost', 'synthesis_yield')

all_inputs = (chitosan_cost, polysorbate_20_cost, azodicarboxamide_97_cost, methylpentane2_cost, 
              glacial_acetic_acid_cost, glutaraldehyde_50_cost, hotplate_cost, hotplate_energy, 
              NaOH_cost, DI_water_cost, chitosan_amount, glacial_acetic_acid_amount, process_temperature,
              NaOH_amount, overnight_settling_time, DI_water_syn, glutaraldehyde_amount, DI_water_rinsing, 
              methylpentane2_amount, tween20_PCPM_amount, azocarboxamide_amount, tween20_PCPT, chitosan_heat_time, 
              electricity_cost, synthesis_yield)

dfinputs = pd.DataFrame({k:v.flatten() for k,v in zip(all_inputs_names, all_inputs)})

df_spearmans_PCPM_cost = (dfinputs.corrwith(df_cost.PCPM_cost, method='spearman')).abs()
df_spearmans_PCPA_cost = (dfinputs.corrwith(df_cost.PCPA_cost, method='spearman')).abs()
df_spearmans_PCPT_cost = (dfinputs.corrwith(df_cost.PCPT_cost, method='spearman')).abs()

df_spearmans_AP_PCPM_impacts = (dfinputs.corrwith(df_impact_AP.PCPM_impacts, method='spearman')).abs()
df_spearmans_AP_PCPA_impacts = (dfinputs.corrwith(df_impact_AP.PCPA_impacts, method='spearman')).abs()
df_spearmans_AP_PCPT_impacts = (dfinputs.corrwith(df_impact_AP.PCPT_impacts, method='spearman')).abs()

df_spearmans_GWP_PCPM_impacts = (dfinputs.corrwith(df_impact_GWP.PCPM_impacts, method='spearman')).abs()
df_spearmans_GWP_PCPA_impacts = (dfinputs.corrwith(df_impact_GWP.PCPA_impacts, method='spearman')).abs()
df_spearmans_GWP_PCPT_impacts = (dfinputs.corrwith(df_impact_GWP.PCPT_impacts, method='spearman')).abs()

df_spearmans_EFW_PCPM_impacts = (dfinputs.corrwith(df_impact_EFW.PCPM_impacts, method='spearman')).abs()
df_spearmans_EFW_PCPA_impacts = (dfinputs.corrwith(df_impact_EFW.PCPA_impacts, method='spearman')).abs()
df_spearmans_EFW_PCPT_impacts = (dfinputs.corrwith(df_impact_EFW.PCPT_impacts, method='spearman')).abs()

df_spearmans_EP_PCPM_impacts = (dfinputs.corrwith(df_impact_EP.PCPM_impacts, method='spearman')).abs()
df_spearmans_EP_PCPA_impacts = (dfinputs.corrwith(df_impact_EP.PCPA_impacts, method='spearman')).abs()
df_spearmans_EP_PCPT_impacts = (dfinputs.corrwith(df_impact_EP.PCPT_impacts, method='spearman')).abs()

df_spearmans_HTC_PCPM_impacts = (dfinputs.corrwith(df_impact_HTC.PCPM_impacts, method='spearman')).abs()
df_spearmans_HTC_PCPA_impacts = (dfinputs.corrwith(df_impact_HTC.PCPA_impacts, method='spearman')).abs()
df_spearmans_HTC_PCPT_impacts = (dfinputs.corrwith(df_impact_HTC.PCPT_impacts, method='spearman')).abs()

df_spearmans_HTNC_PCPM_impacts = (dfinputs.corrwith(df_impact_HTNC.PCPM_impacts, method='spearman')).abs()
df_spearmans_HTNC_PCPA_impacts = (dfinputs.corrwith(df_impact_HTNC.PCPA_impacts, method='spearman')).abs()
df_spearmans_HTNC_PCPT_impacts = (dfinputs.corrwith(df_impact_HTNC.PCPT_impacts, method='spearman')).abs()

df_spearmans_ODP_PCPM_impacts = (dfinputs.corrwith(df_impact_ODP.PCPM_impacts, method='spearman')).abs()
df_spearmans_ODP_PCPA_impacts = (dfinputs.corrwith(df_impact_ODP.PCPA_impacts, method='spearman')).abs()
df_spearmans_ODP_PCPT_impacts = (dfinputs.corrwith(df_impact_ODP.PCPT_impacts, method='spearman')).abs()

df_spearmans_PMPF_PCPM_impacts = (dfinputs.corrwith(df_impact_PMPF.PCPM_impacts, method='spearman')).abs()
df_spearmans_PMPF_PCPA_impacts = (dfinputs.corrwith(df_impact_PMPF.PCPA_impacts, method='spearman')).abs()
df_spearmans_PMPF_PCPT_impacts = (dfinputs.corrwith(df_impact_PMPF.PCPT_impacts, method='spearman')).abs()

df_spearmans_MIR_PCPM_impacts = (dfinputs.corrwith(df_impact_MIR.PCPM_impacts, method='spearman')).abs()
df_spearmans_MIR_PCPA_impacts = (dfinputs.corrwith(df_impact_MIR.PCPA_impacts, method='spearman')).abs()
df_spearmans_MIR_PCPT_impacts = (dfinputs.corrwith(df_impact_MIR.PCPT_impacts, method='spearman')).abs()

dfinputs.to_excel(writer, sheet_name= 'inputs')

df_spearmans_PCPM_cost.to_excel(writer, sheet_name= 'spearmans_PCPM_cost')
df_spearmans_PCPA_cost.to_excel(writer, sheet_name= 'spearmans_PCPA_cost')
df_spearmans_PCPT_cost.to_excel(writer, sheet_name= 'spearmans_PCPT_cost')


df_spearmans_AP_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_AP')
df_spearmans_AP_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_AP')
df_spearmans_AP_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_AP')


df_spearmans_GWP_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_GWP')
df_spearmans_GWP_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_GWP')
df_spearmans_GWP_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_GWP')

df_spearmans_EFW_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_EFW')
df_spearmans_EFW_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_EFW')
df_spearmans_EFW_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_EFW')

df_spearmans_EP_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_EP')
df_spearmans_EP_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_EP')
df_spearmans_EP_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_EP')

df_spearmans_HTC_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_HTC')
df_spearmans_HTC_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_HTC')
df_spearmans_HTC_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_HTC')

df_spearmans_HTNC_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_HTNC')
df_spearmans_HTNC_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_HTNC')
df_spearmans_HTNC_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_HTNC')

df_spearmans_ODP_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_ODP')
df_spearmans_ODP_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_ODP')
df_spearmans_ODP_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_ODP')

df_spearmans_PMPF_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_PMPF')
df_spearmans_PMPF_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_PMPF')
df_spearmans_PMPF_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_PMPF')

df_spearmans_MIR_PCPM_impacts.to_excel(writer, sheet_name= 'spearmans_PCPM_MIR')
df_spearmans_MIR_PCPA_impacts.to_excel(writer, sheet_name= 'spearmans_PCPA_MIR')
df_spearmans_MIR_PCPT_impacts.to_excel(writer, sheet_name= 'spearmans_PCPT_MIR')

spearmans_results = pd.DataFrame(dict(
    MPCP_cost = df_spearmans_PCPM_cost, APCP_cost = df_spearmans_PCPA_cost, TPCP_cost = df_spearmans_PCPT_cost,
    MPCP_AP = df_spearmans_AP_PCPM_impacts , APCP_AP = df_spearmans_AP_PCPA_impacts, TPCP_AP = df_spearmans_AP_PCPT_impacts ,
    MPCP_GWP = df_spearmans_GWP_PCPM_impacts, APCP_GWP = df_spearmans_GWP_PCPA_impacts, TPCP_GWP = df_spearmans_GWP_PCPT_impacts,
    MPCP_EFW = df_spearmans_EFW_PCPM_impacts, APCP_EFW = df_spearmans_EFW_PCPA_impacts, TPCP_EFW = df_spearmans_EFW_PCPT_impacts,
    MPCP_EP = df_spearmans_EP_PCPM_impacts, APCP_EP = df_spearmans_EP_PCPA_impacts, TPCP_EP = df_spearmans_EP_PCPT_impacts,
    MPCP_HTC = df_spearmans_HTC_PCPM_impacts, APCP_HTC = df_spearmans_HTC_PCPA_impacts, TPCP_HTC = df_spearmans_HTC_PCPT_impacts,
    MPCP_HTNC = df_spearmans_HTNC_PCPM_impacts, APCP_HTNC = df_spearmans_HTNC_PCPA_impacts, TPCP_HTNC = df_spearmans_HTNC_PCPT_impacts,
    MPCP_ODP = df_spearmans_ODP_PCPM_impacts, APCP_ODP = df_spearmans_ODP_PCPA_impacts, TPCP_ODP = df_spearmans_ODP_PCPT_impacts,
    MPCP_PMPF = df_spearmans_PMPF_PCPM_impacts, APCP_PMPF = df_spearmans_PMPF_PCPA_impacts, TPCP_PMPF = df_spearmans_PMPF_PCPT_impacts,
    MPCP_MIR = df_spearmans_MIR_PCPM_impacts, APCP_MIR = df_spearmans_MIR_PCPA_impacts, TPCP_MIR = df_spearmans_MIR_PCPT_impacts,
    )).reset_index()

#
spearmans_results.to_excel(writer, sheet_name= 'all_spearmans')

X = list(spearmans_results)
X.remove('index')

Y = list(spearmans_results['index'])

spearmans_results_scatter = pd.melt(spearmans_results, id_vars=['index'], value_vars=X)

#%% ouput data to excel file
writer.save()



#%% plots


#set colors
colors = ['lightgreen', 'orange', 
          'lightblue']


#cost
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
cost_plot = ax.boxplot(df_cost, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(cost_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in cost_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Cost \n [USD/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan Beads", fontname="Arial", fontsize=12)

plt.tick_params(direction='inout')
plt.savefig("cost.png", dpi=1000, bbox_inches="tight")



#AP
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
AP_plot = ax.boxplot(df_impact_AP, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(AP_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in AP_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Acidification potential \n [kg SO2 eq/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')

plt.savefig("impact_AP.png", dpi=1000, bbox_inches="tight")

#GWP
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
GWP_plot = ax.boxplot(df_impact_GWP, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(GWP_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in GWP_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Global warming potential \n [kg CO2 eq/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')

plt.savefig("impact_GWP.png", dpi=1000, bbox_inches="tight")


#EFW
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
EFW_plot = ax.boxplot(df_impact_EFW, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(EFW_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in EFW_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Ecotoxicity: freshwater \n [CTUe/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')

plt.savefig("impact_EFW.png", dpi=1000, bbox_inches="tight")



#EP
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
EP_plot = ax.boxplot(df_impact_EP, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(EP_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in EP_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Eutrophication potential \n [kg N eq/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')

plt.savefig("impact_EP.png", dpi=1000, bbox_inches="tight")



#HTC
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
HTC_plot = ax.boxplot(df_impact_HTC, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(HTC_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in HTC_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Human toxicity: carcinogenic \n [CTUh/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))


plt.tick_params(direction='inout')
plt.savefig("impact_HTC.png", dpi=1000, bbox_inches="tight")



#HTNC
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
HTNC_plot = ax.boxplot(df_impact_HTNC, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(HTNC_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in HTNC_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Human toxicity: non-carcinogenic \n [CTUh/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')

plt.savefig("impact_HTNC.png", dpi=1000, bbox_inches="tight")



#ODP
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
ODP_plot = ax.boxplot(df_impact_ODP, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(ODP_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in ODP_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Ozone depletion potential \n [kg CFC-11 eq/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')
plt.savefig("impact_ODP.png", dpi=1000, bbox_inches="tight")


#PMPF
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
PMPF_plot = ax.boxplot(df_impact_PMPF, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(PMPF_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in PMPF_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Particulate matter formation \n [kg PM2.5 eq/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

plt.tick_params(direction='inout')

plt.savefig("impact_PMPF.png", dpi=1000, bbox_inches="tight")


#MIR
fig = plt.figure(figsize =(3, 3))
ax = fig.add_subplot(111)

# Creating axes instance
MIR_plot = ax.boxplot(df_impact_MIR, patch_artist = True, widths=0.6, flierprops={'markersize': 2})


for patch, color in zip(MIR_plot['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['MPCP', 'APCP', 'TPCP']) 

for median in MIR_plot['medians']:
    median.set(color ='black',linewidth = 1)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Change the y axis label to Arial
ax.set_ylabel("Maximum incremental reactivity \n [kg O3 eq/g]", fontname="Arial", fontsize=12) 

# Change the x axis label to Arial
ax.set_xlabel("Citosan beads", fontname="Arial", fontsize=12)

plt.tick_params(direction='inout')
#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))



plt.savefig("impact_MIR.png", dpi=1000, bbox_inches="tight")




# spearmans
spearmans_results_scatter['value'] = np.where(spearmans_results_scatter['value']<0.05, np.nan, spearmans_results_scatter['value'])
spearmans_results_scatter["value"] *= 300

spearmans_results_scatter['syn_methods'] = spearmans_results_scatter['variable'].str[:4]

strings_to_remove=('electricity_cost', 'DI_water_syn', 'overnight_settling_time', 'process_temperature', 'glacial_acetic_acid_amount', 'NaOH_cost', 'hotplate_cost', 'glacial_acetic_acid_cost', 'methylpentane2_cost', 'azodicarboxamide_97_cost')

spearmans_results_scatter_filtered = spearmans_results_scatter[~spearmans_results_scatter['index'].isin(strings_to_remove)]

def convert_to_color(variable_name):
    color_mapping = dict(zip(spearmans_results_scatter['syn_methods'].unique(), colors))
    return color_mapping.get(variable_name, 'gray')  # Default to gray for unknown variables

# Add a new column 'Color' to the DataFrame
spearmans_results_scatter['color']= spearmans_results_scatter['syn_methods'].apply(convert_to_color)


#data filtering
strings_to_remove=('electricity_cost', 'DI_water_syn', 'overnight_settling_time', 'process_temperature', 'glacial_acetic_acid_amount', 'NaOH_cost', 'hotplate_cost', 'glacial_acetic_acid_cost', 'methylpentane2_cost', 'azodicarboxamide_97_cost')

spearmans_results_scatter_filtered = spearmans_results_scatter[~spearmans_results_scatter['index'].isin(strings_to_remove)]


fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
plt.xticks(rotation=90)

# Set the font name for axis tick labels to be Comic Sans
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

ax.grid(color = 'gray', linestyle = '--', linewidth = 0.5)

ax.set_xlabel("Sustainability indicators", fontname="Arial", fontsize=12)
ax.set_ylabel("Assumptions", fontname="Arial", fontsize=12) 
ax.scatter('variable','index', s='value', alpha=0.8, c='color', data=spearmans_results_scatter_filtered)
plt.tick_params(direction='inout')
plt.margins(y=.04)
plt.savefig("spearmans.png", dpi=1000, bbox_inches="tight")
plt.savefig("spearmans.eps", dpi=1000, bbox_inches="tight")


























