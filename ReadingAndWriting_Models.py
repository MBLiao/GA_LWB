from pathlib import Path
from cobra.io import load_json_model, save_json_model, load_matlab_model, \
    save_matlab_model, read_sbml_model, write_sbml_model, load_model
import logging

data_dir = Path(r"C:\Users\廖文彬\.conda\envs\COBRAenv\Lib\site-packages\cobra\data")
# print(data_dir)
data_dir = data_dir.resolve()

print("mini test files: ")
print(", ".join(str(i) for i in data_dir.glob('mini*.*')))

textbook_model = load_model("textbook")
ecoli_model = load_model("iJO1366")
logging.getLogger("cobra.io.sbml").setLevel(logging.ERROR)
salmonella_model = load_model("salmonella")


# mini_fbc2_path = data_dir / "mini_fbc2.xml"
# read_sbml_model(str(mini_fbc2_path.resolve()))

mini_json_path = data_dir / "mini.json"
load_json_model(str(mini_json_path.resolve()))
save_json_model(textbook_model, "test.json")












