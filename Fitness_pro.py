import pandas as pd
from pathlib import Path
import numpy as np
import warnings
import time
import os
import sys
import matplotlib.pyplot as plt
from ModelOperation import *
import ray
from joblib import Parallel, delayed

warnings.filterwarnings("ignore")
# start_time = time.time()

#
# def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
#     # 这里放置计算个体适应度的代码，假设返回适应度值
#     fitness = fitness_fun(individual=individual, model0=model0, target_table=target_table,
#                           GeneProteinMap=GeneProteinMap, tartgetRxn=tartgetRxn, ref=ref, metrxn=metrxn)
#     return fitness
#
#
# def parallel_fitness_evaluation(population, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
#     # 使用joblib的Parallel和delayed来并行处理
#     fitness_scores = Parallel(n_jobs=-1)\
#         (delayed(evaluate_individual)(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn)
#          for individual in population)
#     return fitness_scores


def fitness_fun(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        for i, gene in enumerate(individual):
            if gene != 0:
                target_gene = target_table.iloc[i]["gene_name"]
                target_enzyme = GeneProteinMap[target_gene]
                # 根据操作调整模型
                if gene == 1:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.000000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.000000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                elif gene == 2:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                elif gene == 3:
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
        try:

            ft1 = time.time()
            solution_m = MPMA.moma(model=model, solution=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantYield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['r_2111']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['r_1714_REV']
                print(growth, product, glucose, mutantYield_m)
                print('-----------------------------------------------------------------------------------------------')

                # calculate the product and biomass fold change in mutant strain compared with wild strain

            else:
                mutantYield_m = 0.01
        except:
            mutantYield_m = 0.01

    return mutantYield_m















