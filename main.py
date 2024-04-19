import pandas as pd
import ModelOperation
from pathlib import Path
from cobra.io import load_matlab_model
from cobra.flux_analysis import flux_variability_analysis
import numpy as np
import warnings
import time
import os
import sys
import matplotlib.pyplot as plt
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

warnings.filterwarnings("ignore")
# start_time = time.time()


def loadModel():
    data_dir = Path("D:\models")
    data_dir = data_dir.resolve()

    ecYeastGEM_path = data_dir / "ecYeastGEM_batch.mat"
    model = load_matlab_model(str(ecYeastGEM_path.resolve()))
    model.solver = 'gurobi'
    model.solver.configuration.verbosity = 0

    return model


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


def WildModel(model):
    # here we firstly simulate the real conditions

    # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

    model.reactions.get_by_id('r_1589').bounds = (2.5, 2.5)  # set the initial rate of the production
    model.objective = 'r_2111'
    solution0 = model.optimize()
    growth0 = solution0.objective_value
    print('growth0: ', growth0)
    print(model.reactions.get_by_id("r_1714_REV").bounds)
    print('glu------0：', solution0.fluxes['r_1714_REV'])
    # fix growth0 and minimization the protein pool
    model.reactions.get_by_id('r_2111').bounds = (0.90 * growth0, growth0)

    model.objective = {model.reactions.prot_pool_exchange: -1}
    solutionWild = model.optimize()

    print('*' * 12)
    print('WILD: ', solutionWild.fluxes['r_1589'], solutionWild.fluxes['r_1714_REV'])
    wildYield = solutionWild.fluxes['r_1589'] / solutionWild.fluxes['r_1714_REV']
    print('wildYield: ', wildYield)
    print('*' * 12)

    # model.objective = 'r_1589'
    # result = model.optimize()
    # print('glu------1：', result.fluxes['r_1714_REV'])

    print('-'*120)

    return model


def PredictionOE(model0, tartgetRxn, changeRxn, OE=4):
    print('-'*120)
    # glucose is unlimited
    model = model0.copy()
    # model.objective = {model.reactions.prot_pool_exchange: -1}
    basic_solution = model.optimize()
    wildYield = basic_solution.fluxes['r_1589'] / basic_solution.fluxes['r_1714_REV']
    enzymeUsage = basic_solution.fluxes[changeRxn]
    print(changeRxn)
    print('wildYield: ', wildYield)
    print('enzymeUsage: ', enzymeUsage)

    enzyme_range = list(frange(1, OE, 0.1))  # for the over-expression


    # model.objective = 'r_2111'
    #
    # growth = []
    # for i in enzyme_range:
    #     model.reactions.get_by_id(changeRxn).lower_bound = enzymeUsage * i
    #     try:
    #         solution = model.optimize()
    #         growth.append(solution.objective_value)
    #     except:
    #         growth.append(None)
    # growth = [x for x in growth if x is not None]
    # growth_min = min(growth)
    # print('MIN growth: ', growth_min)
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
        print('max_target: ', solution2.objective_value)
        model.reactions.get_by_id(tartgetRxn).bounds = (0.999 * solution2.objective_value, solution2.objective_value)
        # we put it in the mutant model
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solutionMut = model.optimize()

        print('.' * 12)
        print('MUT: ', solutionMut.fluxes[tartgetRxn], solutionMut.fluxes['r_1714_REV'])
        print('.' * 12)

        # calculate the product and biomass fold change in mutant strain compared with wild strain
        mutantYield = solutionMut.fluxes[tartgetRxn] / solutionMut.fluxes['r_1714_REV']
    except:
        print('Too high enzyme usage!!')
        mutantYield = None

    print('mutantYield: ', mutantYield)
    print('@')
    print('-'*120)

    return wildYield, mutantYield


def fitness_fun(individual, model, target_table, G, O):
    # 克隆一份模型以进行修改
    modified_model = model.copy()

    # 遍历基因编码
    for i, gene in enumerate(individual):
        if gene == 1:
            # 获取目标基因/酶的名称和调整方式
            target_gene = target_table.iloc[i]["gene_name"]
            operation = target_table.iloc[i]["operation"]

            related_reactions = []
            for reaction in modified_model.reactions:
                for ge in reaction.genes:
                    if str(target_gene) == str(ge):
                        related_reactions.append(reaction)
                        break

            print("Number of related reactions: ", len(related_reactions), "PO.", i)

            # 修改这些反应的上界
            for reaction in related_reactions:
                # 根据操作调整模型
                if operation == "OE":
                    # 上调反应上界
                    # reaction.upper_bound *= 100
                    reaction.lower_bound *= 10
                    # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

                elif operation == "KD":
                    # 下调反应上界
                    reaction.upper_bound *= 0.5
                    # print("########################################")

                elif operation == "KO":
                    reaction.upper_bound = 0
                    reaction.lower_bound = 0
                    # print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")

    # 运行FBA
    T1 = time.time()
    solution = modified_model.optimize()
    T2 = time.time()
    print('@@@@@@@', T2 -T1)
    # 假设适应度是生长速率
    # 检查求解器状态
    if solution.status == 'optimal':
        # 如果优化成功，获取目标值
        fitness = 100 * (solution.fluxes['r_1589']/O) * (solution.fluxes['r_2111']/G)
    else:
        # 如果优化失败，可以设置一个默认值或处理错误
        print("模型优化失败，状态:", solution.status)
        fitness = 0
    print("适应度:", fitness)

    return fitness

def saveExcel(infile, outfile):
    '''
    function to save the dataframe into xlsx format
    :param infile:
    :param outfile:
    :return:
    '''
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer, 'Sheet1')
    writer.save()




ecYeast = loadModel()

rxnID = []
for i, x in enumerate(ecYeast.reactions):
    # print(x)
    rxnID.append(x.id)
rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]
print(len(rxnID_protein_draw))

geneAll = []
for gene in ecYeast.genes:
    # print(gene.id)
    geneAll.append(gene.id)
# get the gene list for the related proreins with kcat number
gene_with_Kcat = [x.replace('draw_prot_', '') for x in rxnID_protein_draw]

# rxnID_protein_draw = ['draw_prot_P61829', 'draw_prot_P38858', 'draw_prot_P32327']

WildModel = WildModel(ecYeast)

# loops for the over-expression
all_yield_w = []
all_yield_m = []
test_protein_draw = []
i = 1
for j in rxnID_protein_draw:
    print('choose protein need over expression: ' + str(j) + '--------------------------')
    test_protein_draw.append(j)
    s = PredictionOE(model0=WildModel, tartgetRxn='r_1589', changeRxn=j)
    all_yield_w.append(s[0])
    all_yield_m.append(s[1])
    i += 1
    if i == 201:
        break

# put the result into a dataframe
overExpression_result = pd.DataFrame({'gene': test_protein_draw, 'yield_wild': all_yield_w,
                                      'yield_mutant': all_yield_m})

saveExcel(overExpression_result, r"C:\Users\Liao\Desktop\ecYeast_OE_output.xls")

print('1')









