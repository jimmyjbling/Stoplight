from rdkit import Chem
from .pk_calculator import CLASSIFICATION_DICT
from .assay_liability_calculator import OUTCOME_DICT

CLASSIFICATION_DICT.update(OUTCOME_DICT)

# inverts the inner dicts of the classification dicts so outcome can be converted to numbers
INVERSE_CLASS_DICT = {key2: {val: key for key, val in val2.items()} for key2, val2 in CLASSIFICATION_DICT.items()}

MOLECULE_PROPERTY_SCORES = {
    'Solubility in Water (mg/L)': ([50, 10], "bigger"),
    'LogP': ([2, 3], "smaller"),
    'Molecular Weight': ([400, 500], "smaller"),
    'Polar Surface Area': ([120, 140], "smaller"),
    'Number of Rotatable Bonds': ([7, 11], "smaller"),
    'FSP3': ([0.3, 0.2], "bigger"),
    "Putative aggregator": 1,
    "Non-aggregator": 0,
    "Possible Interference": 1,
    "No Interference": 0
}

COLORS = {
    0: "green",
    1: "#b09b2a",
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


def get_stoplight(smiles, options):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    results = {}
    overall_scores = []
    # ESOL
    if "Solubility in Water (mg/L)" in options.keys():
        from .esol_calculator import calc_esol
        res = calc_esol(smiles, options)
        res = {key: val + [get_prop_score(key, val[1])] for key, val in res.items()}
        overall_scores += [_[-1] for _ in res.values()]
        res = {key: val + [COLORS[val[-1]]] for key, val in res.items()}
        results.update(res)

    # RDKit Properties
    if any([_ in options.keys() for _ in ['LogP','Molecular Weight','Polar Surface Area','Number of Rotatable Bonds','FSP3']]):
        from .prop_calculator import get_props
        res = get_props(smiles, options)
        res = {key: val + [get_prop_score(key, val[1])] for key, val in res.items()}
        overall_scores += [_[-1] for _ in res.values()]
        res = {key: val + [COLORS[val[-1]]] for key, val in res.items()}
        results.update(res)

    # Get assay liability
    if any([_ in options.keys() for _ in ['Firefly Luciferase interference','Nano Luciferase interference','Redox interference','Thiol interference','AmpC β-lactamase aggregation','Cysteine protease cruzain aggregation']]):
        from .assay_liability_calculator import main as assay_main
        res = assay_main(smiles, **options)
        res = {key: val + [MOLECULE_PROPERTY_SCORES[val[1]]] for key, val in res.items()}
        overall_scores += [_[-1] for _ in res.values()]
        res = {key: val + [COLORS[val[-1]]] for key, val in res.items()}
        results.update(res)

    # Get pk properties
    if any([_ in options.keys() for _ in ['Hepatic Stability', 'Microsomal Half-life Sub-cellular', 'Microsomal Half-life Tissue', 'Renal Clearance', 'BBB Permeability', 'CNS Activity', 'CACO2', 'Plasma Protein Binding', 'Plasma Half-life', 'Microsomal Intrinsic Clearance', 'Oral Bioavailability']]):
        from .pk_calculator import main as pk_main
        results.update(pk_main(smiles, **options))

    if len(overall_scores) > 0:
        overall_score = sum(overall_scores) / len(overall_scores)
    else:
        overall_score = None

    return results, overall_score
