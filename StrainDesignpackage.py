import cobra
import straindesign as sd
import pandas as pd
from pathlib import Path
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

sys.path.append('/path/to/ch/javasoft')

warnings.filterwarnings("ignore")
model = cobra.io.load_model('e_coli_core')

data_dir = Path("D:\models")
data_dir = data_dir.resolve()

ecYeastGEM_path = data_dir / "ecYeastGEM_batch.mat"
model = load_matlab_model(str(ecYeastGEM_path.resolve()))
model.solver = 'gurobi'
model.solver.configuration.verbosity = 0

solution = sd.fba(model)

print(f"Maximum growth: {solution.objective_value}.")
solution = sd.fba(model, obj='draw_prot_P11154', obj_sense='minimize')

print(f"Minimum flux through draw_prot_P11154: {solution.objective_value}.")









