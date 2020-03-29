#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

# Takes the data in GSE110137_rpkms.gct units of reads per kilobase per million mapped reads,
# deletes the column 'NAME', and adds a title to each axis
frame = pd.read_csv('../data/GSE110137_rpkms.gct', sep='\t', skiprows=2, index_col=1)
del frame['NAME']
frame.columns.names, frame.index.names = ['Antibiotics'], ['Genes']

# Genes that are unnamed are removed
frame = frame[np.all([frame.index != '-', frame.index != 'not available'], axis=0)]
frame


# In[2]:


# Sample IDs under !series_sample_id in the GSE110137_series_matrix.txt file
cols = '"PDD_P2_18"	"PDD_P2_44"	"PDD_P2_70"	"PDD_P2_20"	"PDD_P2_46"	"PDD_P2_72"	"PDD_P2_11"	"PDD_P2_37"	"PDD_P2_63"	"PDD_P2_05"	"PDD_P2_31"	"PDD_P2_57"	"PDD_P2_21"	"PDD_P2_47"	"PDD_P2_73"	"PDD_P2_07"	"PDD_P2_33"	"PDD_P2_59"	"PDD_P2_15"	"PDD_P2_41"	"PDD_P2_67"	"PDD_P2_01"	"PDD_P2_02"	"PDD_P2_27"	"PDD_P2_28"	"PDD_P2_53"	"PDD_P2_54"	"PDD_P2_09"	"PDD_P2_35"	"PDD_P2_61"	"PDD_P2_23"	"PDD_P2_49"	"PDD_P2_75"	"PDD_P2_10"	"PDD_P2_36"	"PDD_P2_62"	"PDD_P2_14"	"PDD_P2_40"	"PDD_P2_66"	"PDD_P2_04"	"PDD_P2_30"	"PDD_P2_56"	"PDD_P2_03"	"PDD_P2_29"	"PDD_P2_55"	"PDD_P2_06"	"PDD_P2_32"	"PDD_P2_58"	"PDD_P2_24"	"PDD_P2_50"	"PDD_P2_76"	"PDD_P2_12"	"PDD_P2_38"	"PDD_P2_64"	"PDD_P2_26"	"PDD_P2_52"	"PDD_P2_78"'
cols = cols.replace('"', '').split('\t')

# Compounds tested under !Sample_title from the GSE110137_series_matrix.txt file
compounds = '"AZ-LolCDE"	"AZ-LolCDE"	"AZ-LolCDE"	"B-01"	"B-01"	"B-01"	"Ceftriaxone"	"Ceftriaxone"	"Ceftriaxone"	"Chloramphenicol"	"Chloramphenicol"	"Chloramphenicol"	"Ciprofloxacin"	"Ciprofloxacin"	"Ciprofloxacin"	"Clarithromycin"	"Clarithromycin"	"Clarithromycin"	"Colistin"	"Colistin"	"Colistin"	"DMSO"	"DMSO"	"DMSO"	"DMSO"	"DMSO"	"DMSO"	"Doxycycline"	"Doxycycline"	"Doxycycline"	"Globomycin"	"Globomycin"	"Globomycin"	"Levofloxacin"	"Levofloxacin"	"Levofloxacin"	"Mecillinam"	"Mecillinam"	"Mecillinam"	"Meropenem"	"Meropenem"	"Meropenem"	"Nitrofurantoin"	"Nitrofurantoin"	"Nitrofurantoin"	"Nitroxolin"	"Nitroxolin"	"Nitroxolin"	"Norfloxacin"	"Norfloxacin"	"Norfloxacin"	"PolymyxineB"	"PolymyxineB"	"PolymyxineB"	"Trimethoprim"	"Trimethoprim"	"Trimethoprim"'
compounds = compounds.replace('"', '').split('\t')

# Renames the titles in the original DataFrame with the compound names
col_dict = dict(zip(cols, compounds))
rpkm = frame.rename(columns=col_dict)
rpkm


# In[3]:


from scipy.stats import zscore

# Calcuates the z-score for each gene relative to all other genes in a sample
z = rpkm.apply(zscore).T
z


# In[4]:


import matplotlib

# Sets font as Arial
matplotlib.rcParams['font.sans-serif'] = 'Arial'


# In[5]:


import matplotlib.pyplot as plt
import seaborn as sns

# Clusters each condition by euclidean distance. Please see
# https://seaborn.pydata.org/generated/seaborn.clustermap.html
# for more documentation
h = sns.clustermap(z, figsize=(13, 10), col_cluster=False, cmap='RdBu_r', yticklabels=True, xticklabels=False,
                   cbar_kws={'pad': 0.19, 'shrink': 0.7, 'label': 'Z-score'},
                   cbar_pos=None, cbar=True, vmin=-0.5, vmax=0.5)

# Removes all tick lines in the figure
for i in h.fig.axes:
    i.tick_params(size=0)

# Saves figure to as clustermap_test.png
plt.savefig('../figures/clustermap_test.png', dpi=300)
plt.show()

