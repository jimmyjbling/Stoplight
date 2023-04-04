import argparse

import pandas as pd

from Stoplight.main import get_csv_from_smiles

CONVERT_OPTIONS = {
    'fsp3': 'FSP3',
    'logp': 'LogP',
    'mw': 'Molecular Weight',
    'rot_bonds': 'Number of Rotatable Bonds',
    'psa': 'Polar Surface Area',
    'esol': 'Solubility in Water (mg/L)',
    'blac_agg': 'AmpC β-lactamase aggregation',
    'cprot_agg': 'Cysteine protease cruzain aggregation',
    'fluc_inter': 'Firefly Luciferase interference',
    'nluc_inter': 'Nano Luciferase interference',
    'redox_inter': 'Redox interference',
    'thiol_inter': 'Thiol interference',
    'bbb': 'BBB Permeability',
    'caco2': 'CACO2',
    'cns': 'CNS Activity',
    'hep_stab': 'Hepatic Stability',
    'micro_hf_sub': 'Microsomal Half-life Sub-cellular',
    'micro_hf_t': 'Microsomal Half-life Tissue',
    'micro_clr': 'Microsomal Intrinsic Clearance',
    'o_avail': 'Oral Bioavailability',
    'pla_hf': 'Plasma Half-life',
    'pla_pb': 'Plasma Protein Binding',
    'ren_clr': 'Renal Clearance'
}

PROPERTY_LITERALS = {
    "all": ['fsp3', 'logp', 'mw', 'rot_bonds', 'psa', 'esol', 'blac_agg', 'cprot_agg', 'fluc_inter', 'nluc_inter',
            'redox_inter', 'thiol_inter', 'bbb', 'caco2', 'cns', 'hep_stab', 'micro_hf_sub', 'micro_hf_t', 'micro_clr',
            'o_avail', 'pla_hf', 'pla_pb', 'ren_clr'],
    "molecular": ['fsp3', 'logp', 'mw', 'rot_bonds', 'psa', 'esol'],
    "pk": ['bbb', 'caco2', 'cns', 'hep_stab', 'micro_hf_sub', 'micro_hf_t', 'micro_clr',
           'o_avail', 'pla_hf', 'pla_pb', 'ren_clr'],
    "liability": ['blac_agg', 'cprot_agg', 'fluc_inter', 'nluc_inter', 'redox_inter', 'thiol_inter'],
    'fsp3': ['fsp3'],
    'logp': ['logp'],
    'mw': ['mw'],
    'rot_bonds': ['rot_bonds'],
    'psa': ['psa'],
    'esol': ['esol'],
    'blac_agg': ['blac_agg'],
    'cprot_agg': ['cprot_agg'],
    'fluc_inter': ['fluc_inter'],
    'nluc_inter': ['nluc_inter'],
    'redox_inter': ['redox_inter'],
    'thiol_inter': ['thiol_inter'],
    'bbb': ['bbb'],
    'caco2': ['caco2'],
    'cns': ['cns'],
    'hep_stab': ['hep_stab'],
    'micro_hf_sub': ['micro_hf_sub'],
    'micro_hf_t': ['micro_hf_t'],
    'micro_clr': ['micro_clr'],
    'o_avail': ['o_avail'],
    'pla_hf': ['pla_hf'],
    'pla_pb': ['pla_pb'],
    'ren_clr': ['ren_clr']
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--infile", type=str, required=True,
                        help="path to input file of SMILES")
    parser.add_argument("--smi_col", required=False, default=0,
                        help="index/header of column in csv containing SMILES. Use index if no header (starting at 0)")
    parser.add_argument("--props", nargs="*", required=False, default=("all",),
                        help="properties to calculated for scoring. Names of properties that can passed be found in "
                             "'stoplight_properties_help.txt'. Besides individual property names can pass 'all' for "
                             "all, 'molecular' for molecular properties, 'pk' for pharmacokinetic predictions and "
                             "'liability' for assay liability predictions")
    parser.add_argument("--outfile", type=str, required=False, default="stoplight_results.csv",
                        help="path to write output file to")
    parser.add_argument("--drop_invalid", action="store_true", default=False, required=False,
                        help="drop invalid SMILES from file. Defaults to keeping invalid and listing columns as "
                             "'invalid SMILES'")

    args = parser.parse_args()

    options = {key: True for key in list(set([__ for _ in args.props for __ in PROPERTY_LITERALS[_]]))}
    options["drop_invalid"] = args.drop_invalid

    try:
        _col = int(args.smi_col)
        smiles = pd.read_csv(args.infile, header=None)[0]
    except ValueError:
        smiles = pd.read_csv(args.infile)[args.smi_col]

    csv_text = get_csv_from_smiles(smiles_list=smiles.to_list(), options=options)
