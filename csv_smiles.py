from rdkit.Chem import MolFromSmiles
from io import StringIO
import csv
from Stoplight.main import get_stoplight, INVERSE_CLASS_DICT


def get_csv_from_smiles(smiles_list, options):
    # CSV writer expects a file object, not a string.
    # StringIO can be used to store a string as a file-like object.

    options["make_prop_img"] = False  # do not need to create images for csv

    headers = [key for key, val in options.items() if ((key not in ["make_prop_img", "precision"]) and val)]

    # Add prop headers for datapoints with AD
    prob_headers = [f"{key}_prob" for key in headers if key not in ['LogP', 'Molecular Weight', 'Polar Surface Area', 'Number of Rotatable Bonds', 'FSP3', 'Solubility in Water (mg/L)']]
    # Add AD headers for datapoints with AD
    ad_headers = [f"{key}_AD" for key in headers if key not in ['LogP', 'Molecular Weight', 'Polar Surface Area', 'Number of Rotatable Bonds', 'FSP3', 'Solubility in Water (mg/L)']]

    headers = headers + prob_headers + ad_headers

    string_file = StringIO()
    writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers, "OverallScore"])
    writer.writeheader()

    # loop through all the smiles
    for smiles in smiles_list:
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        # reject invalid smiles without crashing
        if molecule is None:
            row['SMILES'] = f"(invalid){smiles}"
            writer.writerow(row)
            continue

        props, overall_score = get_stoplight(smiles, options)

        props = [[key] + val for key, val in props.items()]

        # loop through all properties for the smiles
        for prop_name, _, pred, pred_proba, ad, score, _ in props:
            if prop_name not in ['LogP', 'Molecular Weight', 'Polar Surface Area', 'Number of Rotatable Bonds', 'FSP3', 'Solubility in Water (mg/L)']:
                print(prop_name, _, pred, pred_proba, ad, score, _)
                row[prop_name] = INVERSE_CLASS_DICT[prop_name][pred]
                try:
                    row[prop_name + "_prob"] = float(pred_proba[:-1]) / 100  # covert back to 0-1 float
                except ValueError:
                    row[prop_name] = "NA"  # if pred_proba is string skip
                    row[prop_name + "_prob"] = "NA"  # covert back to 0-1 float
                row[prop_name + "_AD"] = ad
            else:
                row[prop_name] = pred

        row["OverallScore"] = overall_score
        writer.writerow(row)

    return string_file.getvalue()
