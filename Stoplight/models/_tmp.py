from xgboost import XGBClassifier
from rdkit.Chem.PandasTools import LoadSDF
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

ENDPOINT_DICT = {
    "agg_betalac.sdf": {
        'Non-Aggregator': 0,
        'Putative Aggregator': 1
    },
    "agg_cruzain.sdf": {
        'Non-Aggregator': 0,
        'Putative Aggregator': 1
    },
    "Firefly_Luciferase_counter_assay_training_set_curated.sdf": {
        '0': 0,
        '1': 1
    },
    "Nano_Luciferase_counter_assay_training_set_curated.sdf": {
        '0': 0,
        '1': 1
    },
    "Redox_training_set_curated.sdf": {
        '0': 0,
        '1': 1
    },
    "Thiol_training_set_curated.sdf": {
        '0': 0,
        '1': 1
    }
}

import os
import pickle

for file in os.listdir("C:\\Users\\welln\\OneDrive\\TropshaLab\\Website\\Stoplight\\Stoplight\\models"):
    if not file.endswith(".sdf"):
        continue
    if file in ["agg_betalac.sdf", "agg_cruzain.sdf"]:
        continue
    print("Loading file")
    df = LoadSDF(file)
    print("Making fingerprints")
    df["SMILES"] = df["ROMol"].apply(Chem.MolToSmiles)
    df["fp"] = df["ROMol"].apply(AllChem.GetMorganFingerprintAsBitVect, radius=3, nBits=2048)

    print("Flattening data")
    fp = np.array([list(x) for x in df["fp"]])
    y = np.array([ENDPOINT_DICT[file][x] for x in df["Outcome"]])

    print("Fitting Model")
    model = XGBClassifier(n_jobs=-1)
    model.fit(fp, y)

    print(f"{file} Done")
    with open(f"{file.split('.')[0]}_model.pkl", "bw") as f:
        pickle.dump(model, f)
