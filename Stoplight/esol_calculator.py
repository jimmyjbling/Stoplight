from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from collections import namedtuple
import numpy as np
import pickle
import os

from Stoplight.constants import MODULE_DIR

print(os.getcwd())
MODEL = pickle.load(open(os.path.join(MODULE_DIR, "models", "esol_model.pkl"), 'rb'))
AROMATIC_QUERY = Chem.MolFromSmarts("a")
Descriptor = namedtuple("Descriptor", "mw logp rotors ap")


def calc_ap(mol):
    """
    Calculate aromatic proportion #aromatic atoms/#atoms total
    :param mol: input molecule
    :return: aromatic proportion
    """
    matches = mol.GetSubstructMatches(AROMATIC_QUERY)
    return len(matches) / mol.GetNumAtoms()


def calc_esol_descriptors(mol):
    """
    Calcuate mw,logp,rotors and aromatic proportion (ap)
    :param mol: input molecule
    :return: named tuple with descriptor values
    """
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rotors = Lipinski.NumRotatableBonds(mol)
    ap = calc_ap(mol)
    return Descriptor(mw=mw, logp=logp, rotors=rotors, ap=ap)


def calc_esol_delaney(mol):
    """
    Calculate ESOL based on descriptors in the Delaney paper, coefficients refit for the RDKit using the
    routine refit_esol below
    :param mol: input molecule
    :return: predicted solubility
    """
    intercept = 0.26121066137801696
    coef = {'mw': -0.0066138847738667125, 'logp': -0.7416739523408995, 'rotors': 0.003451545565957996, 'ap': -0.42624840441316975}
    desc = calc_esol_descriptors(mol)
    esol = intercept + coef["logp"] * desc.logp + coef["mw"] * desc.mw + coef["rotors"] * desc.rotors \
            + coef["ap"] * desc.ap

    mwt = Descriptors.MolWt(mol) * 1000

    esol = 10**esol * mwt

    return esol


def calc_esol_xgbmodel(mol):
    desc_funcs = [x[1] for x in Descriptors.descList if x[0] != "Ipc"]
    X = np.array([func(mol) for func in desc_funcs]).reshape(1, -1)

    try:
        esol = MODEL.predict(X)
        esol = float(esol[0])
        mwt = Descriptors.MolWt(mol) * 1000
        esol = (10 ** esol) * mwt

    except Exception:
        esol = calc_esol_delaney(mol)

    return esol


def calc_esol(smiles, options):
    def get_value(value_func, mol):
        return round(value_func(mol), int(options['precision']))

    mol = Chem.MolFromSmiles(smiles)
    return {'Solubility in Water (mg/L)': [0, get_value(calc_esol_xgbmodel, mol), "", ""]}
