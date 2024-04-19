from pathlib import Path
from cobra.io import load_matlab_model
import numpy as np
import warnings
import time
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

warnings.filterwarnings("ignore")
start_time = time.time()
data_dir = Path("D:\models")
data_dir = data_dir.resolve()

ecYeastGEM_path = data_dir / "ecYeastGEM.mat"
model = load_matlab_model(str(ecYeastGEM_path.resolve()))

end_time_1 = time.time()

solution = model.optimize()
print(solution)
print('objective expression', model.objective.expression)
print('objective direction', model.objective.direction)
print('complete model: ', solution.objective_value)
print('-'*120)


for reaction in model.reactions[:5]:
    with model as model:
        reaction.knock_out()
        model.optimize()
        print('%s blocked (bounds: %s), new growth rate %f' %
              (reaction.id, str(reaction.bounds), model.objective.value))
print('-'*120)

end_time_2 = time.time()

model.solver = 'gurobi'
# 创建一个列表来存储不同基因敲除组合的结果
results = []

# 循环十次，每次敲除五个基因
for i in range(10):
    # 从模型中选择要敲除的随机的5个基因
    genes_to_knock_out = np.random.choice(model.genes, 5, replace=False)

    # 敲除选定的基因
    with model:
        for gene in genes_to_knock_out:
            gene.knock_out()
        solution = model.optimize()
        if solution.objective_value is not None:
            results.append((genes_to_knock_out, solution.objective_value))

    # 输出每次敲除的基因组合和对应结果值
    print(
        f"Iteration {i + 1}: Genes Knocked Out : {[gene.id for gene in genes_to_knock_out]},"
        f" Objective Value : {solution.objective_value}")

# 找到最佳结果
best_result = max(results, key=lambda x: x[1])
# 输出最佳结果
print("\n最佳基因敲除组合:", [gene.id for gene in best_result[0]])
print("最佳优化结果:", best_result[1])
print(len(model.genes))
print('-'*120)

end_time_3 = time.time()


# # 创建一个列表来存储单基因敲除结果
# deletion_results = []
#
# # 循环遍历每个基因
# with model:
#     for gene in model.genes:
#         gene.knock_out()
#         solution = model.optimize()
#
#         # 检查目标值是否为None，如果不是，则将结果添加到列表
#         if solution.objective_value is not None:
#             deletion_results.append((gene, solution.objective_value))
#
# # 检查是否有有效的结果
# if deletion_results:
#     # 找到目标值最大的单基因敲除方案
#     best_result_1 = max(deletion_results, key=lambda x: x[1])
#     # 输出最佳结果
#     print("\n最佳单基因敲除:", best_result_1[0].id)
#     print("最佳优化结果:", best_result_1[1])
# else:
#     print("未找到有效的单基因敲除方案")
# print(deletion_results)
# print('-'*120)

# deletion_results = single_gene_deletion(model)
# print(deletion_results)

# start = time.time()  # start timer()
# double_gene_deletion(
#     model, model.genes[:25], processes=2)
# t1 = time.time() - start
# print("Double gene deletions for 200 genes completed in "
#       "%.2f sec with 2 cores" % t1)
#
# start = time.time()  # start timer()
# double_gene_deletion(
#     model, model.genes[:25], processes=1)
# t2 = time.time() - start
# print("Double gene deletions for 200 genes completed in "
#       "%.2f sec with 1 core" % t2)
#
# print("Speedup of %.2fx" % (t2 / t1))

end_time_4 = time.time()

print(f"读取模型时间：{end_time_1 - start_time} s")
print(f"简单优化时间：{end_time_2 - end_time_1} s")
print(f"组合基因敲除时间：{end_time_3 - end_time_2} s")
print(f"单基因敲除遍历时间：{end_time_4 - end_time_3} s")
print('-'*120)




Re = model.reactions.get_by_id('r_1264')
Re_1 = model.reactions.get_by_id('r_1265')

print(Re.reaction)
print(Re.lower_bound, "< REACTION_succinate <", Re.upper_bound)
print(Re.reversibility)

print(Re_1.reaction)
print(Re_1.lower_bound, "< REACTION_succinate <", Re_1.upper_bound)
print(Re_1.reversibility)
print('-'*120)

Suc = model.metabolites.get_by_id('s_1458')
Suc_1 = model.metabolites.get_by_id('s_1459')
Suc_2 = model.metabolites.get_by_id('s_1460')
print(Suc.name)
print(Suc.formula)
print(Suc_1.name)
print(Suc_1.formula)
print(Suc_2.name)
print(Suc_2.formula)
print('-'*120)

print(model.summary())
print(model.metabolites.s_1458.summary())
print('-'*120)
print(model.metabolites.s_1459.summary())
print('-'*120)
print(model.metabolites.s_1460.summary())
print('-'*120)


from cobra.util.solver import linear_reaction_coefficients

with model:
    model.objective = 'r_1265'
    print('objective in this context:', model.objective.expression)
    print('objective direction', model.objective.direction)

    linear_reaction_coefficients(model)
    print(model.solver.status())
    print(model.optimize().objective_value)
    print('-' * 120)

