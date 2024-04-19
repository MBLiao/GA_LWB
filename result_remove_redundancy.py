import pandas as pd
from pathlib import Path
import numpy as np
import warnings
import time
import os
import sys
import matplotlib.pyplot as plt
from ModelOperation import *
from Fshow import *

warnings.filterwarnings("ignore")


pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 1000)
np.set_printoptions(linewidth=400, threshold=200)  # np.inf表示正无穷
warnings.filterwarnings("ignore")

ecYeast = loadModel()

# wildSolution_moma = WildModel_MOMA(ecYeast)
wildSolution_moma = WildModel_Growth(ecYeast)
# wildSolution_fba = WildModel_FBA(ecYeast)

ecYeast.reactions.get_by_id('r_2111').lower_bound = 0.1
ecYeast.reactions.get_by_id('r_1589').lower_bound = 0.1
ecYeast.reactions.get_by_id('r_1714_REV').lower_bound = 0.1

print(ecYeast.reactions.get_by_id('r_2111').bounds)
print(ecYeast.reactions.get_by_id('r_1589').bounds)
print(ecYeast.reactions.get_by_id('r_1714_REV').bounds)

wildYield = wildSolution_moma.fluxes['r_1589'] / wildSolution_moma.fluxes['r_1714_REV']

rxnID = []
for i, x in enumerate(ecYeast.reactions):
    # print(x)
    rxnID.append(x.id)
rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]
geneAll = []
for gene in ecYeast.genes:
    # print(gene.id)
    geneAll.append(gene.id)

metabolic_solution = pd.Series()
for rid in wildSolution_moma.fluxes.index:
    if rid.startswith('r_'):
        metabolic_solution[rid] = wildSolution_moma.fluxes[rid]

enzyme_solution = pd.Series()
for rid in wildSolution_moma.fluxes.index:
    if rid.startswith('draw_'):
        enzyme_solution[rid] = wildSolution_moma.fluxes[rid]

arm_solution = pd.Series()
for rid in wildSolution_moma.fluxes.index:
    if rid.startswith('arm_'):
        arm_solution[rid] = wildSolution_moma.fluxes[rid]

# get the gene list for the related proteins with kcat number
gene_with_Kcat = [x.replace('draw_prot_', '') for x in rxnID_protein_draw]
print('----------------------------MODEL INFORMATION----------------------------')
print('rxnID: ', len(rxnID))
print('geneAll: ', len(geneAll))
print('rxnID_protein_draw: ', len(rxnID_protein_draw))
print('gene_with_Kcat: ', len(gene_with_Kcat))
print('metabolic reactions: ', len(metabolic_solution))
print('enzyme reactions: ', len(enzyme_solution))
print('arm reactions: ', len(arm_solution))
print('WILD YEILD: ', wildYield)
print(rxnID_protein_draw[:5])
print(gene_with_Kcat[:5])
print('-------------------------------------------------------------------------')

GeneProteinMap = {}
df = pd.read_excel('/Users/liaowenbin/Documents/models/ecModel_batch.xls')
GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

ecFSEOF_tabel = pd.read_csv('/Users/liaowenbin/Desktop/ecFSEOF_MAPPING.csv')
# iBridge_tabel = pd.read_csv(r"C:\Users\Liao\Desktop\ibridge_2PE_result.csv")
# ecFactory_tabel = pd.read_csv(r"C:\Users\Liao\Desktop\ecYeast_2PE_ecFactory_output.csv")

ind = [2, 0, 1, 2, 0, 0, 0, 2, 0, 2, 3, 3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 1,
       0, 1, 3, 0, 3, 3, 0, 0, 2, 1, 2, 3, 2, 0, 1, 0, 2, 3, 1, 2, 0, 0, 0, 2, 1, 3]
# ind = [0, 1, 0, 0, 3, 0, 1, 0, 1, 1, 0, 2, 0, 0, 3, 0, 3, 0, 0, 0, 2, 0, 2, 0, 0, 0,
#        0, 0, 0, 2, 2, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 2, 0, 0, 1, 0, 0, 2, 0]
# ind = [2, 0, 1, 2, 0, 0, 0, 2, 0, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ind = [1, 1, 3, 2, 0, 0, 0, 3, 0, 1, 2, 1, 2, 1, 1, 1, 0, 2, 3, 0, 3, 3, 0, 0, 0, 3,
       1, 1, 2, 1, 1, 3, 0, 2, 2, 2, 2, 0, 0, 0, 0, 2, 2, 1, 2, 3, 0, 2, 2, 3, 0, 0]
ind = [0, 1, 0, 0, 3, 0, 1, 0, 1, 1, 0, 2, 0, 0, 3, 0, 3, 0, 0, 0, 2, 0, 2, 0, 0, 0,
       0, 0, 0, 2, 2, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 2, 0, 0, 1, 0, 0, 2, 0]
ind = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0]
ans = fs_fun(ind, ecYeast, ecFSEOF_tabel, GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)
print(ans)
for p in range(len(ind)):
    with ind as d:
        if d[p] != 0:
            d[p] = 0
            be = fs_fun(ind, ecYeast, ecFSEOF_tabel, GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)
            if be != ans:
                print(p)





























