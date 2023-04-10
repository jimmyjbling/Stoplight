from flask import Flask, render_template, request, abort, Response, jsonify

from Stoplight.main import is_smiles
from smiles import get_molecule_data_from_smiles
from csv_smiles import get_csv_from_smiles
from Stoplight.assay_liability_calculator import MODEL_DICT as AL_MODEL_DICT
from Stoplight.prop_calculator import MOLECULE_PROPERTIES
from Stoplight.pk_calculator import MODEL_DICT as PK_MODEL_DICT


ALL_PROPS = {"mp": sorted(list(MOLECULE_PROPERTIES.keys())) + ["Solubility in Water (mg/L)"],
             "al": sorted(list(AL_MODEL_DICT.values())),
             "pk": sorted(list(PK_MODEL_DICT.keys()))}

app = Flask(__name__)


@app.route('/')
def home():
    print("main")
    return render_template('index.html')


@app.route('/score-cutoff')
def score_cutoff():
    print("cutoff")
    return render_template('score_table.html')


@app.route('/score-description')
def score_description():
    print("cutoff")
    return render_template('score_cutoffs.html')


@app.route('/endpoints', methods=['GET'])
def mol_properties():
    print("endpoints")
    return jsonify(ALL_PROPS), 200


@app.route('/check_smiles', methods=['POST'])
def check_smiles():
    _tmp = {"is_smile": int(is_smiles(request.json.get('smiles')))}
    print(_tmp)
    return _tmp


@app.route('/smiles', methods=['POST'])
def smiles():
    print("smiles")

    data = get_molecule_data_from_smiles(request.json.get('smiles'), request.json.get('options'))

    if data is None:
        return abort(400)

    return data


@app.route('/smiles-csv', methods=['POST'])
def smiles_csv():
    print("csv")

    _smiles = request.json.get('smiles')

    if len(_smiles) > 1000:
        return abort(413)

    print("here")
    csv = get_csv_from_smiles(_smiles, request.json.get('options'))

    return Response(
        csv,
        mimetype="text/csv",
        headers={"Content-disposition": "attachment;"}
    )

