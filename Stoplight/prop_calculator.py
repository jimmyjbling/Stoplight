from collections import OrderedDict

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED


MOLECULE_PROPERTIES = OrderedDict({
    'LogP': Crippen.MolLogP,
    'Molecular Weight': Descriptors.MolWt,
    'Polar Surface Area': lambda mol: QED.properties(mol).PSA,
    'Number of Rotatable Bonds': Lipinski.NumRotatableBonds,
    'FSP3': Lipinski.FractionCSP3
})


def get_props(smiles, options):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    def get_value(value_func):
        return round(value_func(mol), int(options['precision']))

    mol_props = {prop_name: [0, get_value(value_func), "", ""] for (prop_name, value_func) in MOLECULE_PROPERTIES.items() if prop_name in options}

    return mol_props
