import cobra
import straindesign as sd
import pandas as pd
from pathlib import Path
from cobra.io import load_matlab_model
from cobra.flux_analysis import moma
import numpy as np
import warnings
import time
import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

data_dir = Path("/Users/liaowenbin/Documents/models")
data_dir = data_dir.resolve()
ecYeastGEM_path = data_dir / "ecYeastGEM_batch.mat"
model0 = load_matlab_model(str(ecYeastGEM_path.resolve()))
model0.solver = 'gurobi'
model0.solver.configuration.verbosity = 0
cobra.Configuration().tolerance = 1e-9

with model0 as model:
    rr = model.reactions.get_by_id('r_2111').upper_bound
    for i, x in enumerate(model.reactions):
        if model.reactions.get_by_id(x.id).upper_bound == rr:
            model.reactions.get_by_id(x.id).upper_bound = 1000
        # print(model.reactions.get_by_id(x.id).bounds)
    print(f"COBRApy\'s solver is \'{cobra.Configuration().solver.__name__}',"
          f" StrainDesign selects {sd.select_solver()}.")
    print('---------------------------------')

    solution1 = sd.fba(model, constraints='r_1589>=2.0')
    print('Growth: ', solution1.fluxes['r_2111'])
    print('Product: ', solution1.fluxes['r_1589'])
    print('Glucose: ', solution1.fluxes['r_1714_REV'])
    print('Enzyme: ', solution1.fluxes['prot_pool_exchange'])
    print('---------------------------------')


    solution2 = sd.fba(model, obj='r_1714_REV', obj_sense='minimize',
                       constraints=[f'r_2111>={solution1.objective_value*0.9}', 'r_1589>=2.0'])
    print('Growth: ', solution2.fluxes['r_2111'])
    print('Product: ', solution2.fluxes['r_1589'])
    print('Glucose: ', solution2.fluxes['r_1714_REV'])
    print('Enzyme: ', solution2.fluxes['prot_pool_exchange'])
    print('---------------------------------')

    solution3 = sd.fba(model, obj='prot_pool_exchange', obj_sense='minimize',
                       constraints=[f'r_2111>={solution1.objective_value*0.9}', 'r_1589>=2.0',
                                    f'r_1714_REV<={solution2.objective_value*1.1}'])

    print('Growth: ', solution3.fluxes['r_2111'])
    print('Product: ', solution3.fluxes['r_1589'])
    print('Glucose: ', solution3.fluxes['r_1714_REV'])
    print('Enzyme: ', solution3.fluxes['prot_pool_exchange'])
    print('---------------------------------')

    # fba_sol = sd.fba(model)
    # pfba1_sol = sd.fba(model, solver='glpk', pfba=1)
    # print(f"The sum of fluxes of the regular FBA: {sum([abs(v) for v in fba_sol.fluxes.values()])} " + \
    #       f"is usually higher than of the parsimoniuos FBA: {sum([abs(v) for v in pfba1_sol.fluxes.values()])}")

    # df = pd.DataFrame(list(pfba1_sol.fluxes.items()), columns=['Reaction', 'Flux Value'])
    # print(pfba1_sol.fluxes.values())
    # print(df.describe())
    # print(f"Maximum growth FBA: {fba_sol.objective_value}.")
    # print(f"Maximum growth pFBA: {pfba1_sol.objective_value}.")
    #
    # cons_fva = [f'r_2111>={solution1.objective_value*0.9}',
    #             'r_1589>=2.0',
    #             f'r_1714_REV<={solution2.objective_value*1.1}']
    # cons_fva = [f'r_2111>=0.1',
    #             'r_1589>=0.1']
    # cons_fva = []
    # solution_fva = sd.fva(model, constraints=cons_fva, processes=1)
    # print(solution_fva)
    # solution_fva.to_csv(r"C:\Users\Liao\Desktop\1000FVA_0.9_2_1.1.csv", index=True)


    # numerator = 'r_1589'
    # denominator = 'r_1714_REV'
    # constraint = ['r_2111 >= 0.1', 'r_1589 >= 2', 'r_1714_REV >= 0.1']
    #
    # sol_fba_g = sd.fba(model, obj='r_2111', obj_sense='maximize', constraints=constraint)
    # sol_fba_p = sd.fba(model, obj='r_1589', obj_sense='maximize', constraints=constraint)
    # fba_yield_g = sol_fba_g.fluxes[numerator] / sol_fba_g.fluxes['r_1714_REV']
    # fba_yield_p = sol_fba_p.fluxes[numerator] / sol_fba_p.fluxes['r_1714_REV']
    #
    # sol_yopt = sd.yopt(model, obj_num=numerator, obj_den=denominator, obj_sense='maximize', constraints=constraint)
    # yopt_yield = sol_yopt.objective_value
    #
    # print(f"(Maximum growth) yield (FBA): {fba_yield_g}.")
    # print(f"(Maximum product) yield (FBA): {fba_yield_p}.")
    # print(f"Maximum yield (yOpt): {yopt_yield}.")


    GeneProteinMap = {}
    df = pd.read_excel('/Users/liaowenbin/Documents/models/ecModel_batch.xls')
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

    ecFSEOF_tabel = pd.read_csv('/Users/liaowenbin/Desktop/ecFSEOF_MAPPING.csv')
    individual = [2, 0, 1, 2, 0, 0, 0, 2, 0, 2, 3, 3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 1,
                  0, 1, 3, 0, 3, 3, 0, 0, 2, 1, 2, 3, 2, 0, 1, 0, 2, 3, 1, 2, 0, 0, 0, 2, 1, 3]
    individual = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0]
    target_table = ecFSEOF_tabel
    ref = solution1

    with model as MT:
        # 遍历基因编码
        for i, gene in enumerate(individual):
            if gene != 0:
                # 获取目标基因/酶的名称和调整方式
                target_gene = target_table.iloc[i]["gene_name"]

                target_enzyme = GeneProteinMap[target_gene]
                # print(target_gene, target_enzyme)
                # 根据操作调整模型
                if gene == 1:
                    print('E' * 12)
                    print('E' * 12)
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.000000001:
                        MT.reactions.get_by_id(target_enzyme).lower_bound = 0.000000004
                    else:
                        MT.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                    print(MT.reactions.get_by_id(target_enzyme).bounds)
                elif gene == 2:
                    print('D' * 12)
                    print('D' * 12)
                    enzymeUsage = ref.fluxes[target_enzyme]
                    MT.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    print(MT.reactions.get_by_id(target_enzyme).bounds)

                elif gene == 3:
                    print('O' * 12)
                    print('O' * 12)
                    MT.reactions.get_by_id(target_enzyme).upper_bound = 0
                    print(MT.reactions.get_by_id(target_enzyme).bounds)
        cons = []
        sd.plot_flux_space(MT, ('r_2111', 'r_1589'), constraints=cons)




    print('1')
