from rdkit.Chem import MolFromSmiles, Draw
from Stoplight.main import get_stoplight


def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None

    print(options)

    props, overall_score = get_stoplight(smiles_str, options)

    props = [[key] + val for key, val in props.items()]

    assay_conensus = []
    for prop_name, _, pred, pred_proba, ad, score, _ in props:
        if prop_name in ["Firefly Luciferase interference", "Nano Luciferase interference", "Redox interference",
                         "Thiol interference", "AmpC β-lactamase aggregation", "Cysteine protease cruzain aggregation"]:
            assay_conensus.append(float(pred_proba[:-1]) / 100 if pred in ["No Interference", "Non-aggregator"] else 1 - (float(pred_proba[:-1]) / 100))

    if len(assay_conensus) > 0:
        _vals = [float(_) for _ in assay_conensus if assay_conensus != "NA"]
        assay_conensus = sum(_vals) / len(_vals) if len(_vals) > 0 else "NA"
    else:
        assay_conensus = -1

    if overall_score == 0:
        with open("./static/green.svg", "r") as f:
            stoplight = f.read()
    elif overall_score <= 1:
        with open("./static/yellow.svg", "r") as f:
            stoplight = f.read()
    elif overall_score <= 2:
        with open("./static/red.svg", "r") as f:
            stoplight = f.read()
    else:
        stoplight = ""

    svg = Draw.MolsToGridImage([molecule], useSVG=True, molsPerRow=1)

    svg = str(svg).replace("'200px'", "'400px'")

    return {
        'svg': svg,
        'SMILES': smiles_str,
        'molProperties': props,
        'assayCon': assay_conensus,
        'stoplight': stoplight,
        'overall': overall_score
    }
