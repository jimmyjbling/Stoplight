from rdkit import Chem

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


def get_stoplight(smiles, options):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    results = {}
    # ESOL
    if "Solubility in Water (mg/L)" in options.keys():
        from esol_calculator import calc_esol
        res = calc_esol(smiles)
        res = {key: [val] + [get_prop_score(key, val)] for key, val in res.items()}
        res = {key: val + [COLORS[val[1]]] for key, val in res.items()}
        results.update(res)

    # RDKit Properties
    if any([_ in options.keys() for _ in ['LogP','Molecular Weight','Polar Surface Area','Number of Rotatable Bonds','FSP3']]):
        from prop_calculator import get_props
        res = get_props(smiles, options)
        res = {key: [val] + [get_prop_score(key, val)] for key, val in res.items()}
        res = {key: val + [COLORS[val[1]]] for key, val in res.items()}
        results.update(res)

    # Get assay liability
    if any([_ in options.keys() for _ in ['Firefly Luciferase interference','Nano Luciferase interference','Redox interference','Thiol interference','AmpC β-lactamase aggregation','Cysteine protease cruzain aggregation']]):
        from assay_liability_calculator import main as assay_main
        res = assay_main(smiles, **options)
        res = {key: val + [MOLECULE_PROPERTY_SCORES[val[0]]] for key, val in res.items()}
        res = {key: val + [COLORS[val[-1]]] for key, val in res.items()}
        results.update(res)

    # for now lets ignore the pk stuff
    # # Get pk properties
    # if any([_ in options.keys() for _ in ['Hepatic Stability', 'Microsomal Half-life Sub-cellular','Microsomal Half-life Tissue', 'Renal Clearance', 'BBB Permeability', 'CNS Activity', 'CACO2', 'Plasma Protein Binding', 'Plasma Half-life', 'Microsomal Intrinsic Clearance', 'Oral Bioavailability']]):
    #     from pk_calculator import main as pk_main
    #     results.update(pk_main(smiles, options))

    return results
