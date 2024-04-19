from cobra.io import load_model
import time

model = load_model("textbook")

solution = model.optimize()
print(solution)
print(solution.objective_value)
print('-'*120)

start_time = time.time()
print(model.optimize().objective_value)
end_time = time.time()
execution_time = end_time - start_time
print(f"代码执行时间：{execution_time} 秒")
print('-'*120)

start_time = time.time()
print(model.slim_optimize())
end_time = time.time()
execution_time = end_time - start_time
print(f"代码执行时间：{execution_time} 秒")
print('-'*120)

print(model.summary())
print('-'*120)

print(model.metabolites.nadh_c.summary())
print('-'*120)

print(model.metabolites.atp_c.summary())
print('-'*120)


'''
Changing the Objectives
'''
biomass_rxn = model.reactions.get_by_id("Biomass_Ecoli_core")


from cobra.util.solver import linear_reaction_coefficients
linear_reaction_coefficients(model)


# change the objective to ATPM
model.objective = "ATPM"
print(model.reactions.get_by_id('ATPM'))
# The upper bound should be 1000, so that we get
# the actual optimal value
model.reactions.get_by_id("ATPM").upper_bound = 1000.
linear_reaction_coefficients(model)
print(model.optimize().objective_value)


'''
FBA will not give always give unique solution, 
because multiple flux states can achieve the same optimum. 
FVA (or flux variability analysis) finds the ranges of each metabolic flux at the optimum
'''
# from cobra.flux_analysis import flux_variability_analysis
# flux_variability_analysis(model, model.reactions[:10])

print('1')





