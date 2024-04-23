import csv
import math
import time
from io import StringIO

from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import MolFromSmiles

from .pk_calculator import CLASSIFICATION_DICT as PK_CLASSIFICATION_DICT
from .assay_liability_calculator import OUTCOME_DICT

CLASSIFICATION_DICT = PK_CLASSIFICATION_DICT.copy()

CLASSIFICATION_DICT.update(OUTCOME_DICT)

ASSAY_LIABILITIES = ["Firefly Luciferase interference", "Nano Luciferase interference", "Redox interference",
                     "Thiol interference", "AmpC β-lactamase aggregation", "Cysteine protease cruzain aggregation"]


# inverts the inner dicts of the classification dicts so outcome can be converted to numbers
INVERSE_CLASS_DICT = {key2: {val: key for key, val in val2.items()} for key2, val2 in CLASSIFICATION_DICT.items()}

MOLECULE_PROPERTY_SCORES = {
    'Solubility in Water (mg/L)': ([50, 10], "bigger"),
    'ALogP': ([2, 3], "smaller"),
    'Molecular Weight': ([400, 500], "smaller"),
    'Polar Surface Area': ([120, 140], "smaller"),
    'Number of Rotatable Bonds': ([7, 11], "smaller"),
    'FSP3': ([0.3, 0.2], "bigger"),
    'Number of Rings': ([2, 1], "bigger"),
    'HBD': ([3, 5], "smaller"),
    'HBA': ([6, 10], "smaller"),
    "Num Heavy Atoms": ([50, 70], "smaller"),
    'Num Saturated Quaternary Carbons': (None, None),
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


def _stoplight_colors(val):
    if val == 0:
        return "green"
    if 0 < val <= 1.0:
        return "yellow"
    if 1.0 < val <= 2.0:
        return "red"


def get_csv_from_smiles(smiles_list, options):
    # CSV writer expects a file object, not a string.
    # StringIO can be used to store a string as a file-like object.

    # print(options)

    options["make_prop_img"] = False  # do not need to create images for csv

    headers = [key for key, val in options.items() if ((key not in ["make_prop_img", "precision"]) and val)]

    headers = headers + [f"{key}_proba" for key in headers if key in PK_CLASSIFICATION_DICT.keys()]

    string_file = StringIO()
    if any([_ in ASSAY_LIABILITIES for _ in headers]):
        writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers, "OverallScore", "StoplightColor"], delimiter=",")
    else:
        writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers, "OverallScore", "StoplightColor"], delimiter=",")
    writer.writeheader()

    # loop through all the smiles
    for smiles in tqdm(smiles_list):
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        # reject invalid smiles without crashing
        if molecule is None:
            if "drop_invalid" not in options.keys() or not options["drop_invalid"]:
                for col in writer.fieldnames:
                    row[col] = "invalid SMILES"
                row['SMILES'] = f"{smiles}"
                writer.writerow(row)
            continue

        props, overall_score = get_stoplight(smiles, options)
        overall_score = round(overall_score, int(options["precision"]))

        props = [[key] + val for key, val in props.items()]
        # loop through all properties for the smiles
        for prop_name, _, pred, pred_proba, ad, score, __ in props:
            if prop_name not in ['ALogP', 'Molecular Weight', 'Polar Surface Area', 'Number of Rotatable Bonds', 'FSP3',
                                 'Solubility in Water (mg/L)', 'Number of Rings', 'HBD', 'HBA', "Num Heavy Atoms",
                                 'Num Saturated Quaternary Carbons']:

                if prop_name in ASSAY_LIABILITIES:
                    print(INVERSE_CLASS_DICT[prop_name][pred])
                    row[prop_name] = round(float(pred_proba[:-1]) / 100.0, int(options["precision"])) if INVERSE_CLASS_DICT[prop_name][pred] == 1 \
                        else round(1 - (float(pred_proba[:-1]) / 100.0), int(options["precision"]))
                else:
                    row[prop_name] = pred
                    try:
                        row[prop_name+"_proba"] = round(float(pred_proba[:-1]) / 100.0, int(options["precision"]))
                    except ValueError:
                        row[prop_name+"_proba"] = pred_proba
                # if prop_name in ["Firefly Luciferase interference", "Nano Luciferase interference", "Redox interference", "Thiol interference", "AmpC β-lactamase aggregation", "Cysteine protease cruzain aggregation"]:
                #     overall_liability_score.append(row[prop_name + "_prob"] if pred in ["No Interference", "Non-aggregator"] else 1-row[prop_name + "_prob"])
            else:
                row[prop_name] = pred

        row["OverallScore"] = overall_score
        row["StoplightColor"] = _stoplight_colors(overall_score)
        # if len(overall_liability_score) > 0:
        #     _vals = [float(_) for _ in overall_liability_score if overall_liability_score != "NA"]
        #     row["ConsensusAssayScore"] = round(1 - (sum(_vals) / len(_vals)), int(options["precision"])) if len(_vals) > 0 else "NA"
        writer.writerow(row)

    return string_file.getvalue()


def get_prop_score(prop_name, prop_val):
    score_range, order = MOLECULE_PROPERTY_SCORES[prop_name]
    if score_range is None:
        return 0
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
    if any([_ in options.keys() for _ in ['LogP', 'Molecular Weight', 'Polar Surface Area', 'Number of Rotatable Bonds', 'FSP3',
                                 'Solubility in Water (mg/L)', 'Number of Rings', 'HBD', 'HBA', "Num Heavy Atoms",
                                 'Num Saturated Quaternary Carbons']]):
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
        overall_score = round(overall_score, int(options["precision"]))

    else:
        overall_score = None

    return results, overall_score


def is_smiles(smi):
    if smi == "":
        return True
    return True if Chem.MolFromSmiles(smi) is not None else False
