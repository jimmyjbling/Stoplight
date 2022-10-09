from rdkit.Chem import MolFromSmiles, Draw
from Stoplight.main import get_stoplight


def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None

    print(options)

    props, overall_score = get_stoplight(smiles_str, options)

    props = [[key] + val for key, val in props.items()]

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
        'stoplight': stoplight,
        'overall': overall_score
    }
