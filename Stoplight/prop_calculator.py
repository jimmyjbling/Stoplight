from collections import OrderedDict

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors
import re

def _num_saturated_quaternary_carbons(mol):
    return len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#6]")))


def _num_consecutive_rings(mol):
    res1 = re.sub(r'[^0-9()]', "", Chem.MolToSmiles(mol))
    res2 = re.sub(r'\(([0-9])\1+\)', "", res1)
    res3 = re.sub(r'[()]', "", res2)

    max_set = 0
    _set = set()
    leader = -1
    for c in res3:
        if leader == -1:
            leader = c
            continue
        if leader == c:
            max_set = max([len(_set), max_set])
            _set = set()
            leader = -1
        else:
            _set.add(c)

    return max_set


MOLECULE_PROPERTIES = OrderedDict({
    'ALogP': Crippen.MolLogP,
    'Molecular Weight': Descriptors.MolWt,
    'Polar Surface Area': lambda mol: QED.properties(mol).PSA,
    'Number of Rotatable Bonds': Lipinski.NumRotatableBonds,
    'FSP3': Lipinski.FractionCSP3,
    'Number of Rings': rdMolDescriptors.CalcNumRings,
    'HBD': rdMolDescriptors.CalcNumHBD,
    'HBA': rdMolDescriptors.CalcNumHBA,
    "Num Heavy Atoms": rdMolDescriptors.CalcNumHeavyAtoms,
    'Num Saturated Quaternary Carbons': _num_saturated_quaternary_carbons
})


def get_props(smiles, options):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    def get_value(value_func):
        return round(value_func(mol), int(options['precision']))

    mol_props = {prop_name: [0, get_value(value_func), "", ""] for (prop_name, value_func) in MOLECULE_PROPERTIES.items() if prop_name in options}

    return mol_props
