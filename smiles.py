from collections import OrderedDict
from rdkit.Chem import MolFromSmiles, Descriptors, Draw, Crippen, Lipinski, QED
from esol_calculator import calc_esol_xgbmodel
from xgboost import XGBRegressor

import pickle
from functools import partial
MODEL = pickle.load(open("./static/esol_model.pkl", 'rb'))

MOLECULE_PROPERTIES = OrderedDict({
    'Predicted Solubility in Water (mg/L)': partial(calc_esol_xgbmodel, model=MODEL),
    'LogP': Crippen.MolLogP,
    'Molecular Weight': Descriptors.MolWt,
    'Polar Surface Area': lambda mol: QED.properties(mol).PSA,
    'Number of Rotatable Bonds': Lipinski.NumRotatableBonds,
    'FSP3': Lipinski.FractionCSP3
})


MOLECULE_PROPERTY_SCORES = OrderedDict({
    'Predicted Solubility in Water (mg/L)': ([50, 10], "bigger"),
    'Pred.ESol': ([50, 10], "bigger"),
    'LogP': ([2, 3], "smaller"),
    'logP': ([2, 3], "smaller"),
    'Molecular Weight': ([400, 500], "smaller"),
    'Mwt': ([400, 500], "smaller"),
    'Polar Surface Area': ([120, 140], "smaller"),
    'PSA': ([120, 140], "smaller"),
    'Number of Rotatable Bonds': ([7, 11], "smaller"),
    'Rot.Bonds': ([7, 11], "smaller"),
    'FSP3': ([0.3, 0.2], "bigger"),
})


COLORS = {
    0: "green",
    1: "gold",
    2: "red"
}


def get_prop_score(prop_name, prop_val):
    score_range, order = MOLECULE_PROPERTY_SCORES[prop_name]
    if order == "smaller":
        score = 2
        if prop_val < score_range[1]:
            score = 1
        if prop_val < score_range[0]:
            score = 0
    else:
        score = 0
        if prop_val < score_range[0]:
            score = 1
        if prop_val < score_range[1]:
            score = 2
    return score


def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None

    print(options)

    def get_value(value_func):
        return round(value_func(molecule), int(options['precision']))

    # Use list of tuples instead of dict to maintain key order
    mol_props = [(prop_name, get_value(value_func)) for (prop_name, value_func) in MOLECULE_PROPERTIES.items()
                 if prop_name in options]

    mol_props = [(prop_name, prop_val, get_prop_score(prop_name, prop_val)) for (prop_name, prop_val) in mol_props]

    mol_props = [(prop_name, prop_val, prop_score, COLORS[prop_score]) for (prop_name, prop_val, prop_score) in mol_props]

    overall_score = sum([x[2] for x in mol_props]) / len(mol_props)

    if overall_score == 0:
        with open("./images/green.svg", "r") as f:
            stoplight = f.read()
    elif overall_score <= 1:
        with open("./images/yellow.svg", "r") as f:
            stoplight = f.read()
    else:
        with open("./images/red.svg", "r") as f:
            stoplight = f.read()

    print(stoplight)

    return {
        'svg': Draw.MolsToGridImage([molecule], useSVG=True, molsPerRow=1),
        'SMILES': smiles_str,
        'molProperties': mol_props,
        'stoplight': stoplight
    }
