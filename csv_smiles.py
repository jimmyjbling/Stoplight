import csv
from rdkit.Chem import MolFromSmiles
from io import StringIO
from smiles import MOLECULE_PROPERTIES, get_prop_score

CONDENSED_NAMES = ['Pred.ESol', 'logP', 'Mwt', 'PSA', 'Rot.Bonds', 'FSP3']


def get_csv_from_smiles(smiles_list, options):
    # CSV writer expects a file object, not a string. 
    # StringIO can be used to store a string as a file-like object.
    string_file = StringIO()
    writer = csv.DictWriter(string_file, fieldnames=['SMILES', *CONDENSED_NAMES, "score", "score breakdown"])
    writer.writeheader()
    for smiles in smiles_list:
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        if molecule is None:
            row['SMILES'] = f"(invalid){smiles}"

        for index, (key, value) in enumerate(MOLECULE_PROPERTIES.items()):
            if key in options:
                row[CONDENSED_NAMES[index]] = round(value(molecule), int(options['precision']))

        row['score'] = sum([get_prop_score(key, val) for key, val in row.items() if key != "SMILES"])
        row['score breakdown'] = ";".join([str(get_prop_score(key, val)) for key, val in row.items() if key not in ["SMILES", "score"]])

        writer.writerow(row)

    return string_file.getvalue()
