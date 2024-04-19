import cobra
from cobra.io import load_model

# "iJO1366" and "salmonella" are also valid arguments
model = load_model("textbook")

print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))
print('-'*120)

print(model)  # a special type of list called a cobra.DictList
print(model.reactions[29])  # get the 30th reaction in the model
for i in range(9):
    print(model.reactions[i])
print(model.metabolites.get_by_id("atp_c"))  # get the cytosolic atp metabolite object (the id is “atp_c”)
print(model.reactions.EX_glc__D_e.bounds)
print('-'*120)


'''
Reactions
'''
pgi = model.reactions.get_by_id("PGI")
print(pgi)
print(pgi.name)
print(pgi.reaction)
print('-'*120)

# reaction upper and lower bounds
print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
print(pgi.reversibility)

old_bounds = pgi.bounds
pgi.bounds = (0, 1000.0)
print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
print("Reversibility after modification:", pgi.reversibility)
pgi.bounds = old_bounds
print("Reversibility after resetting:", pgi.reversibility)
print('-'*120)

old_bounds = pgi.bounds
print('Upper bound prior to setting new lower bound:', pgi.upper_bound)
pgi.upper_bound = 1100
print('Upper bound after setting new upper bound:', pgi.upper_bound)
pgi.bounds = old_bounds
print('-'*120)

# mass balanced
print(pgi.check_mass_balance())

pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
print(pgi.reaction)
print(pgi.check_mass_balance())

pgi.subtract_metabolites({model.metabolites.get_by_id("h_c"): -1})
print(pgi.reaction)
print(pgi.check_mass_balance())

pgi.reaction = "g6p_c --> f6p_c + h_c + green_eggs + ham"
print(pgi.reaction)
pgi.reaction = "g6p_c <=> f6p_c"
print(pgi.reaction)
print('-'*120)


'''
Metabolites
'''
atp = model.metabolites.get_by_id("atp_c")
print(atp)
print(atp.name)
print(atp.compartment)
print(atp.charge)
print(atp.formula)
print(len(atp.reactions))
print(model.metabolites.get_by_id("g6p_c").reactions)
print(model.metabolites.get_by_id("atp_c").reactions)
print('-'*120)


'''
Genes
'''
gpr = pgi.gpr
print(gpr)
gpr_string = pgi.gene_reaction_rule
print(gpr_string)

print(pgi.genes)

pgi_gene = model.genes.get_by_id("b4025")
print(pgi_gene)
print(pgi_gene.reactions)  # Each gene keeps track of the reactions it catalyzes
print('-'*120)

# Altering the gene_reaction_rule will create new gene objects if necessary and update all relationships
pgi.gene_reaction_rule = "(spam or eggs)"
print(pgi.genes)

# Corresponding gene objects also exist. These objects are tracked by the reactions itself, as well as by the model
print(pgi_gene.reactions)
print(model.genes.get_by_id("spam"))
print('-'*120)

# The knock_out_model_genes function will evaluate the GPR
# and set the upper and lower bounds to 0 if the reaction is knocked out.
cobra.manipulation.knock_out_model_genes(model, ["spam"])
print("after 1 KO: %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))

cobra.manipulation.knock_out_model_genes(model, ["eggs"])
print("after 2 KO:  %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))
print('-'*120)

# When knocking out model genes in a context, it is reversed when leaving the context
model = load_model('textbook')
for reaction in model.reactions[:5]:
    with model as model:
        reaction.knock_out()
        model.optimize()
        print('%s blocked (bounds: %s), new growth rate %f' %
              (reaction.id, str(reaction.bounds), model.objective.value))
print([reaction.bounds for reaction in model.reactions[:5]])
print('-'*120)

# Nested contexts are also supported
print('original objective: ', model.objective.expression)
with model:
    model.objective = 'ATPM'
    print('print objective in first context:', model.objective.expression)
    with model:
        model.objective = 'ACALD'
        print('print objective in second context:', model.objective.expression)
    print('objective after exiting second context:', model.objective.expression)
print('back to original objective:', model.objective.expression)

# While it does not have any actual effect, for syntactic convenience it is also
# possible to refer to the model by a different name than outside the context. Such as
# with model as inner:
#     inner.reactions.PFK.knock_out




