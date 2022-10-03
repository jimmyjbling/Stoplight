from flask import Flask, render_template, request, abort, Response, jsonify
from smiles import get_molecule_data_from_smiles, MOLECULE_PROPERTIES
from csv_smiles import get_csv_from_smiles

app = Flask(__name__)


@app.route('/')
def home():
    print("main")
    return render_template('index.html')


@app.route('/score-cutoff')
def score_cutoff():
    print("cutoff")
    return render_template('score_table.html')


@app.route('/endpoints', methods=['GET'])
def mol_properties():
    print("prop")
    return jsonify(list(MOLECULE_PROPERTIES.keys())), 200


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

    csv = get_csv_from_smiles(request.json.get('smiles'), request.json.get('options'))

    return Response(
        csv,
        mimetype="text/csv",
        headers={"Content-disposition": "attachment;"}
    )
