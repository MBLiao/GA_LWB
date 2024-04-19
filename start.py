import numpy as np
import pandas as pd
import ray
from joblib import Parallel, delayed
from GA_loose import *
from Fitness import *
from ModelOperation import *


# 参数设置
POPULATION_SIZE = 100  # 种群大小
GENE_LENGTH = 52      # 基因长度
ACTIVE_GENES = 5       # 活跃基因数
MAX_GENERATIONS = 2000   # 最大迭代次数
CROSSOVER_RATE = 0.8   # 交叉率
MUTATION_RATE = 0.05   # 变异率
# 收敛检测参数
CONVERGENCE_THRESHOLD = 0.01  # 定义适应度变化的阈值
CONVERGENCE_STEPS = 1000        # 连续多少代适应度变化小于阈值时认为已收敛

# Ray初始化
ray.init()


@ray.remote
def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_fun(individual=individual, model0=model0, target_table=target_table,
                          GeneProteinMap=GeneProteinMap, tartgetRxn=tartgetRxn, ref=ref, metrxn=metrxn)
    return fitness


def parallel_fitness_evaluation(population, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    futures = [evaluate_individual.remote(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn)
               for individual in population]
    results = ray.get(futures)
    return results


if __name__ == '__main__':
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
    df = pd.read_excel(r'D:\models\ecModel_batch.xls')
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

    ecFSEOF_tabel = pd.read_csv(r"C:\Users\Liao\Desktop\ecFSEOF_MAPPING.csv")
    # iBridge_tabel = pd.read_csv(r"C:\Users\Liao\Desktop\ibridge_2PE_result.csv")
    # ecFactory_tabel = pd.read_csv(r"C:\Users\Liao\Desktop\ecYeast_2PE_ecFactory_output.csv")

    time_all_start = time.time()

    # 初始化用于保存每代最优适应度和最优个体的变量
    best_fitness_over_generations = []
    best_individual_over_generations = []
    mean_fitness_over_generations = []
    # 初始化收敛检测相关变量
    last_fitness = 0
    convergence_counter = 0

    population = initialize_pop_loose(POPULATION_SIZE, GENE_LENGTH)
    fitness_scores = parallel_fitness_evaluation(population, ecYeast, ecFSEOF_tabel,
                                                 GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)

    # 遗传算法主循环
    for generation in range(MAX_GENERATIONS):
        ti1 = time.time()

        selected_population = roulette_wheel_selection(population, fitness_scores)

        c_pop = []
        c_fitness_scores = []
        for i in range(0, len(selected_population), 2):
            parent1, parent2 = selected_population[i], selected_population[i + 1]
            child1, child2, fg = crossover(parent1, parent2, CROSSOVER_RATE, GENE_LENGTH)
            if fg == 1:
                c_pop.extend([child1, child2])
        if c_pop:
            c_fitness_scores = parallel_fitness_evaluation(c_pop, ecYeast, ecFSEOF_tabel,
                                                           GeneProteinMap, 'r_1589', wildSolution_moma,
                                                           metabolic_solution)
            population.extend(c_pop)
            fitness_scores.extend(c_fitness_scores)

        m_pop = []
        m_fitness_scores = []
        for individual in population:
            mutated_ind, mflag = mutation(individual, MUTATION_RATE, GENE_LENGTH)
            if mflag == 1:
                m_pop.append(mutated_ind)
        if m_pop:
            if len(m_pop) > 1:
                m_fitness_scores = parallel_fitness_evaluation(m_pop, ecYeast, ecFSEOF_tabel,
                                                               GeneProteinMap, 'r_1589', wildSolution_moma,
                                                               metabolic_solution)
            else:
                ms = fitness_fun(m_pop[0], ecYeast, ecFSEOF_tabel,
                                 GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)
                m_fitness_scores.append(ms)
            population.extend(m_pop)
            fitness_scores.extend(m_fitness_scores)

        print('-' * 120)

        print('-' * 120)

        # 找出这一代的最优适应度和最优个体
        max_fitness = max(fitness_scores)
        best_individual = population[fitness_scores.index(max_fitness)]
        mean_fitness = np.mean(fitness_scores)

        # 保存最优适应度和最优个体
        best_fitness_over_generations.append(max_fitness)
        best_individual_over_generations.append(best_individual)
        mean_fitness_over_generations.append(mean_fitness)
        ti2 = time.time()

        # 输出这一代的信息
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print(len(c_pop))
        print(len(m_pop))
        print(len(population))
        print(len(fitness_scores))
        print(f"Generation {generation}: "
              f"Max Fitness = {max_fitness}\n"
              f"Mean Fitness = {mean_fitness}\n"
              f"Time: {ti2 - ti1}\n"
              f"Solvable Num of individual: {sum(i != 0.01 for i in fitness_scores)}")
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        # 收敛检测
        # 如果连续几代的最优适应度没有显著变化，可以提前终止算法
        fitness_change = abs(max_fitness - last_fitness)
        if fitness_change < CONVERGENCE_THRESHOLD:
            convergence_counter += 1
            if convergence_counter >= CONVERGENCE_STEPS:
                print(f"Convergence reached at generation {generation}.")
                break
        else:
            convergence_counter = 0
        last_fitness = max_fitness

        # 更新种群
        population, fitness_scores = select_fixed_size_population(population, fitness_scores, POPULATION_SIZE)

    # 在循环结束后，可以输出或保存整个进程的最优结果
    print(f"Best Fitness over all generations: {max(best_fitness_over_generations)}")
    best_overall = \
        best_individual_over_generations[best_fitness_over_generations.index(max(best_fitness_over_generations))]
    print(f"Best Individual over all generations: {best_overall}")

    data = {'individual': best_individual_over_generations, 'best fitness': best_fitness_over_generations,
            'mean fitness': mean_fitness_over_generations}
    df = pd.DataFrame(data)

    # 将DataFrame保存到CSV文件中
    df.to_csv(r"C:\Users\Liao\Desktop\GA_test_1.csv", index=False)

    # 绘制最大适应度的变化
    plt.plot(best_fitness_over_generations, color='red')

    plt.title('Best Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.show()

    plt.plot(mean_fitness_over_generations, color='blue')
    plt.title('Mean Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.show()

    time_all_end = time.time()
    print('Total Time:', time_all_end - time_all_start)
    ray.shutdown()
    print('\n1')
