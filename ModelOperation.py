import pandas as pd
from pathlib import Path
import cobra
from cobra.io import load_matlab_model
from cobra.flux_analysis import moma
import MPMA
from cobra.flux_analysis import pfba
import numpy as np
import warnings
import time
import os
import sys
import matplotlib.pyplot as plt

import gurobipy as gp


warnings.filterwarnings("ignore")
# start_time = time.time()


def frange(start, stop, step):
    """
    This function is like range, step can be float value, like 0.1
    :param start:
    :param stop:
    :param step:
    :return:
    """
    i = start
    while i < stop:
        yield i
        i += step


def loadModel():
    data_dir = Path("/Users/liaowenbin/Documents/models")
    data_dir = data_dir.resolve()

    ecYeastGEM_path = data_dir / "ecYeastGEM_batch.mat"
    model = load_matlab_model(str(ecYeastGEM_path.resolve()))
    model.solver = 'gurobi'
    model.solver.configuration.verbosity = 0
    cobra.Configuration().tolerance = 1e-9

    return model


def saveExcel(infile, outfile):
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer, 'Sheet1')
    writer.save()


def WildModel_Growth(Model):
    with Model as model:
        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.reactions.get_by_id('r_1589').bounds = (2, 2.5)  # set the initial rate of the production
        model.objective = 'r_2111'
        tg = time.time()
        solution4 = model.optimize()
        tge = time.time()
        growth0 = float(solution4.objective_value)
        print("--------STEP 1--------")
        print('TIME: ', tge - tg)
        print('fluxes[r_2111]: ', solution4.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution4.fluxes['r_1589'])
        print('glucose uptake: ', solution4.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution4.fluxes['prot_pool_exchange'])

        print('growth0: ', growth0)
        print('---------------------------------')

    return solution4


def WildModel_MOMA(Model):
    with Model as model:
        print('--------MODEL INFORMATION--------')
        print('Reactions:', len(model.reactions))
        print('Metabolites:', len(model.metabolites))
        print('Genes', len(model.genes))
        print('Glucose uptake bound: ', model.reactions.get_by_id("r_1714_REV").bounds)
        print('-' * 40)

        # 基本模型求解
        print("--------BASIC RESULT1--------")
        model.objective = 'r_2111'
        solution1 = model.optimize()
        print(solution1)
        print('fluxes[r_2111]: ', solution1.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution1.fluxes['r_1589'])
        print('glucose uptake: ', solution1.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('enzymeUsage(P11154): ', solution1.fluxes['draw_prot_P11154'])
        print('---------------------------------')

        print("--------BASIC RESULT2--------")
        model.objective = 'r_1589'
        solution2 = model.optimize()
        print(solution2)
        print('fluxes[r_2111]: ', solution2.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution2.fluxes['r_1589'])
        print('glucose uptake: ', solution2.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution2.fluxes['prot_pool_exchange'])
        print('enzymeUsage(P11154): ', solution2.fluxes['draw_prot_P11154'])
        print('---------------------------------')

        print("--------BASIC RESULT3--------")
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution3 = model.optimize()
        print(solution3)
        print('fluxes[r_2111]: ', solution3.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution3.fluxes['r_1589'])
        print('glucose uptake: ', solution3.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution3.fluxes['prot_pool_exchange'])
        print('enzymeUsage(P11154): ', solution3.fluxes['draw_prot_P11154'])
        print('---------------------------------')

        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.reactions.get_by_id('r_1589').bounds = (2, 2.5)  # set the initial rate of the production
        model.objective = 'r_2111'
        solution4 = model.optimize()
        growth0 = float(solution4.objective_value)
        print("--------STEP 1--------")
        print('fluxes[r_2111]: ', solution4.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution4.fluxes['r_1589'])
        print('glucose uptake: ', solution4.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution4.fluxes['prot_pool_exchange'])

        print('growth0: ', growth0)
        print('---------------------------------')

        # STEP 2
        # fix growth0 and minimization the glucose uptake
        model.reactions.get_by_id('r_2111').bounds = (0.90 * growth0, growth0)
        model.objective = {model.reactions.get_by_id('r_1714_REV'): -1}
        solution5 = model.optimize()
        glucose0 = -float(solution5.objective_value)
        print("--------STEP 2--------")
        print('fluxes[r_2111]: ', solution5.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution5.fluxes['r_1589'])
        print('glucose uptake: ', solution5.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution5.fluxes['prot_pool_exchange'])

        print('glucose0: ', glucose0)
        print('---------------------------------')

        # STEP 3
        # fix glucose uptake and minimization the protein pool
        model.reactions.get_by_id('r_1714_REV').bounds = (glucose0, 1.1 * glucose0)
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution6 = model.optimize()
        print("--------STEP 3--------")
        print('fluxes[r_2111]: ', solution6.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution6.fluxes['r_1589'])
        print('glucose uptake: ', solution6.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution6.fluxes['prot_pool_exchange'])
        print('---------------------------------')

        print('*' * 40)
        print('WILD: ', solution6.fluxes['r_1589'], solution6.fluxes['r_1714_REV'])
        wildYield = solution6.fluxes['r_1589'] / solution6.fluxes['r_1714_REV']
        BPCY = wildYield * solution6.fluxes['r_2111']
        print('wildYield: ', wildYield)
        print('BPCY_W: ', BPCY)
        print('*' * 40)

        print(model.reactions.get_by_id('r_1714_REV').bounds)
        print(model.reactions.get_by_id('r_2111').bounds)

    return solution6


def WildModel_FBA(Model):
    model = Model.copy()
    with Model as model:
        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.objective = 'r_2111'
        solution4 = model.optimize()
        growth0 = float(solution4.objective_value)
        print("--------STEP 1--------")
        print('fluxes[r_2111]: ', solution4.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution4.fluxes['r_1589'])
        print('glucose uptake: ', solution4.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution4.fluxes['prot_pool_exchange'])

        print('growth0: ', growth0)
        print('---------------------------------')

        # STEP 2
        # fix growth0 and minimization the glucose uptake
        model.reactions.get_by_id('r_2111').bounds = (0.80 * growth0, growth0)
        model.objective = {model.reactions.get_by_id('r_1714_REV'): -1}
        solution5 = model.optimize()
        glucose0 = -float(solution5.objective_value)
        print("--------STEP 2--------")
        print('fluxes[r_2111]: ', solution5.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution5.fluxes['r_1589'])
        print('glucose uptake: ', solution5.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution5.fluxes['prot_pool_exchange'])

        print('glucose0: ', glucose0)
        print('---------------------------------')

        # STEP 3
        # fix glucose uptake and minimization the protein pool
        model.reactions.get_by_id('r_1714_REV').bounds = (glucose0, 1.2 * glucose0)
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution6 = model.optimize()
        protein_pool0 = -float(solution6.objective_value)
        print("--------STEP 3--------")
        print('fluxes[r_2111]: ', solution6.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution6.fluxes['r_1589'])
        print('glucose uptake: ', solution6.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution6.fluxes['prot_pool_exchange'])
        print('---------------------------------')

        # STEP 4
        # fix the protein pool and max product
        model.reactions.get_by_id('prot_pool_exchange').bounds = (protein_pool0, 1.2 * protein_pool0)
        model.objective = {model.reactions.get_by_id('r_1589'): 1}
        solution7 = model.optimize()
        print("--------STEP 4--------")
        print('fluxes[r_2111]: ', solution7.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution7.fluxes['r_1589'])
        print('glucose uptake: ', solution7.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution7.fluxes['prot_pool_exchange'])
        print('---------------------------------')


        print('*' * 40)
        print('WILD: ', solution7.fluxes['r_1589'], solution7.fluxes['r_1714_REV'])
        wildYield = solution7.fluxes['r_1589'] / solution7.fluxes['r_1714_REV']
        BPCY = wildYield * solution7.fluxes['r_2111']
        print('wildYield: ', wildYield)
        print('BPCY_W: ', BPCY)
        print('*' * 40)

        print(model.reactions.get_by_id('r_1714_REV').bounds)
        print(model.reactions.get_by_id('r_2111').bounds)

    return solution7


def KnockOut(model0, tartgetRxn, changeRxn, ref):
    with model0 as model:
        try:
            model.reactions.get_by_id(changeRxn).upper_bound = 0

            solution_m = MPMA.moma(model=model, solution=ref, linear=False)
            distance = solution_m.objective_value
            print('*' * 40)
            print('Distance: ', distance)
            print("--------MOMA--------")
            growth = solution_m.fluxes['r_2111']
            product = solution_m.fluxes[tartgetRxn]
            glucose = solution_m.fluxes['r_1714_REV']
            print(growth, product, glucose)
            print('-' * 20)

            # calculate the product and biomass fold change in mutant strain compared with wild strain
            mutantYield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
            BPCY_m = float(mutantYield_m) * float(solution_m.fluxes['r_2111'])

        except:
            mutantYield_m = 0
            BPCY_m = 0

        print('MUTANTYIELD_m: ', mutantYield_m)
        print('BPCY_m: ', BPCY_m)

        print('*' * 40)

    # return mutantYield_m, mutantYield_p, BPCY_m, BPCY_p, distance
    return mutantYield_m, BPCY_m, distance, growth, product, glucose


def PredictionOE(model0, tartgetRxn, changeRxn, OE=5):

    # glucose is unlimited
    model = model0.copy()
    # here we firstly simulate the real conditions

    # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

    model.reactions.get_by_id(tartgetRxn).bounds = (2.5, 2.5)  # set the initial rate of the production
    model.objective = 'r_2111'
    solution0 = model.optimize()
    growth0 = solution0.objective_value

    # fix growth0 and minimization the protein pool
    model.reactions.get_by_id('r_2111').bounds = (0.90 * growth0, growth0)

    model.objective = {model.reactions.prot_pool_exchange: -1}
    solutionWild = model.optimize()

    print('-'*40)
    print(solutionWild.fluxes[tartgetRxn], solutionWild.fluxes['r_1714_REV'])
    print('-'*40)

    wildYield = solutionWild.fluxes[tartgetRxn] / solutionWild.fluxes['r_1714_REV']
    enzymeUsage = solutionWild.fluxes[changeRxn]
    print('wildYield: ', wildYield)
    print('enzymeUsage: ', enzymeUsage)

    enzyme_range = list(frange(1, OE, 0.1))  # for the over-expression


    # model.objective = 'r_2111'

    # growth = []
    # model.reactions.get_by_id(changeRxn).lower_bound = enzymeUsage * OE
    # try:
    #     solution = model.optimize()
    #     growth = solution.objective_value
    # except:
    #     growth = None
    # growth = [x for x in growth if x is not None]
    # growth_min = min(growth)
    # if growth_min > 0.1:
    #     index = [i for i, x in enumerate(growth) if x == growth_min]
    #     real_factor = enzyme_range[index[0]]
    #     print('real_factor: ', real_factor)

    # mutant model
    # for the mutant model, we over-express the gene and firstly observe whether the max growth is changed
    # over_expression_factor---OE
    model.objective = tartgetRxn
    model.reactions.get_by_id(changeRxn).lower_bound = enzymeUsage * OE * 0.99
    model.reactions.get_by_id("r_2111").lower_bound = 0.1  # fix the minimum growth rate

    # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

    model.reactions.get_by_id(tartgetRxn).bounds = (0, 1000)  # open the target

    try:
        solution2 = model.optimize()
        model.summary()
        print('target: ', solution2.objective_value)
        model.reactions.get_by_id(tartgetRxn).bounds = (0.999 * solution2.objective_value, solution2.objective_value)
        # we put it in the mutant model
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solutionMut = model.optimize()

        print('-' * 120)
        print(solutionMut.fluxes[tartgetRxn], solutionMut.fluxes['r_1714_REV'])
        print('-' * 120)

        # calculate the product and biomass fold change in mutant strain compared with wild strain
        mutantYield = solutionMut.fluxes[tartgetRxn] / solutionMut.fluxes['r_1714_REV']
    except:
        print('Too high enzyme usage!!')
        mutantYield = 0
    print('product fold change:' + str(mutantYield))
    return wildYield, mutantYield


# # loops for the over-expression
# all_yield_w = []
# all_yield_m = []
# all_growth = []
# all_BPCY_m = []
# all_product = []
# all_glucose = []
# all_dis = []
# all_fold = []
# all_status =[]
# test_protein_draw = []
# i = 1
#
# population = ini_pop(POPULATION_SIZE, GENE_LENGTH, ACTIVE_GENES)
#
# for ind in population:
#     print('--Individual Selected: ' + str(ind) + '--')
#     test_protein_draw.append(ind)
#     # s = KnockOut(model0=ecYeast, tartgetRxn='r_1589', changeRxn=j, ref=wildSolution_moma)
#     s = fitness_fun(individual=ind, model0=ecYeast, target_table=ecFSEOF_tabel,
#                     GeneProteinMap=GeneProteinMap, tartgetRxn='r_1589', ref=wildSolution_moma)
#
#     fold = s[0]/wildYield
#     print('--FOLD--: ', fold)
#     print('*' * 40)
#     all_yield_w.append(wildYield)
#     all_yield_m.append(s[0])
#     all_BPCY_m.append(s[1])
#     all_dis.append(s[2])
#     all_growth.append(s[3])
#     all_product.append(s[4])
#     all_glucose.append(s[5])
#     all_status.append(s[6])
#     all_fold.append(fold)
#     # i += 1
#     # if i == 7:
#     #     break
#
# # put the result into a dataframe
# ecFSEOF_result = pd.DataFrame({'gene': test_protein_draw, 'yield_wild': all_yield_w,
#                                 'yield_mutant_m': all_yield_m,
#                                 'growth_mutant': all_growth,
#                                 'product_mutant': all_product,
#                                 'glucose_mutant': all_glucose,
#                                 'BPCY_m': all_BPCY_m,
#                                 'distance': all_dis,
#                                 'fold': all_fold,
#                                 'status': all_status}
#                               )
#
# ecFSEOF_result.to_csv(r"C:\Users\Liao\Desktop\ecFSEOF_realone.csv")








