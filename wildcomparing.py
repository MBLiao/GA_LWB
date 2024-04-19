import numpy as np
import pandas as pd
from ModelOperation import *


if __name__ == '__main__':
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf表示正无穷
    warnings.filterwarnings("ignore")

    ecYeast = loadModel()

    wildSolution_moma = WildModel_MOMA(ecYeast)
    wildSolution_growth = WildModel_FBA(ecYeast)
    # wildSolution_fba = WildModel_FBA(ecYeast)

    ecYeast.reactions.get_by_id('r_2111').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1589').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1714_REV').lower_bound = 0.1

    print(ecYeast.reactions.get_by_id('r_2111').bounds)
    print(ecYeast.reactions.get_by_id('r_1589').bounds)
    print(ecYeast.reactions.get_by_id('r_1714_REV').bounds)

    wildYield = wildSolution_moma.fluxes['r_1589'] / wildSolution_moma.fluxes['r_1714_REV']
    wildYield_growth = wildSolution_growth.fluxes['r_1589'] / wildSolution_growth.fluxes['r_1714_REV']

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
    print('WILD YEILD_growth: ', wildYield_growth)
    print('-------------------------------------------------------------------------')

    print('\n1')
