from cobra import Model, Reaction, Metabolite
model = Model('example_model')

reaction = Reaction('R_3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

ACP_c = Metabolite(
    'ACP_c',
    formula='C11H21N2O7PRS',
    name='acyl-carrier-protein',
    compartment='c')
omrsACP_c = Metabolite(
    'M3omrsACP_c',
    formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment='c')
co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
malACP_c = Metabolite(
    'malACP_c',
    formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein',
    compartment='c')
h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
ddcaACP_c = Metabolite(
    'ddcaACP_c',
    formula='C23H43N2O8PRS',
    name='Dodecanoyl-ACP-n-C120ACP',
    compartment='c')


# SId
# letter   ::=   ’a’..’z’,’A’..’Z’
# digit    ::=   ’0’..’9’
# idChar   ::=   letter | digit | ’_’
# SId      ::=   ( letter | ’_’ ) idChar*

reaction.add_metabolites({
    malACP_c: -1.0,
    h_c: -1.0,
    ddcaACP_c: -1.0,
    co2_c: 1.0,
    ACP_c: 1.0,
    omrsACP_c: 1.0
})

print(reaction.reaction)  # This gives a string representation of the reaction

reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
print(reaction.genes)

# At this point in time, the model is still empty
print(f'{len(model.reactions)} reactions initially')
print(f'{len(model.metabolites)} metabolites initially')
print(f'{len(model.genes)} genes initially')
print('-'*120)

model.add_reactions([reaction])

# The objects have been added to the model
print(f'{len(model.reactions)} reactions')
print(f'{len(model.metabolites)} metabolites')
print(f'{len(model.genes)} genes')
print('-'*120)

# Iterate through the the objects in the model
print("Reactions")
print("---------")
for x in model.reactions:
    print("%s : %s" % (x.id, x.reaction))

print("")
print("Metabolites")
print("-----------")
for x in model.metabolites:
    print('%9s : %s' % (x.id, x.formula))

print("")
print("Genes")
print("-----")
for x in model.genes:
    associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" %
          (x.id, "{" + ", ".join(associated_ids) + "}"))


'''
Last we need to set the objective of the model. Here, 
we just want this to be the maximization of the flux in the 
single reaction we added and we do this by assigning the 
reaction’s identifier to the objective property of the model.
'''
model.objective = 'R_3OAS140'
print(model.objective.expression)
print(model.objective.direction)
print('-'*120)


'''
Model Validation
'''
import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(model, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)

pprint(report)
print('-'*120)


'''
Exchanges, Sinks and Demands
'''
print("exchanges", model.exchanges)
print("demands", model.demands)
print("sinks", model.sinks)

model.add_metabolites([
    Metabolite(
    'glycogen_c',
    name='glycogen',
    compartment='c'
    ),
    Metabolite(
    'co2_e',
    name='CO2',
    compartment='e'
    ),
])

# create exchange reaction
model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")


# boundary reactions
print(model.boundary)
# metabolic reactions
set(model.reactions) - set(model.boundary)
