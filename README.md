# STOPLIGHT

STOPLIGHT Score Calculator: Used to predict/calculate varies properties and assign a hit progression score to a given chemical or list of chemicals. If you use please cite \[INSERT CITATION HERE]. There is a [webserver](https://stoplight.mml.unc.edu/) that runs these models, but for large numbers of compounds, running locally using this code is much more effective

# Requirements
Install the requirements from the requirements.txt file. Additionally, if you want to run the webserver, you need to install flask and gunicorn

# Command line use
After downloading, `Stoplight/stoplight.py` can be called from the command line with `python stoplight.py --help`

`--infile` is required and is the file location for a csv of SMILES to predict properties for. Requires that csv has header and is comma seperated
`--smiles_col` is the name of the column containing the SMILES strings of interest. Defaults to "SMILES"
`--outfile` is the file location of where the output csv file should go. Defaults to `\[current-working-directory\]/stoplight_results.csv`
`--props` is all the properties you want to calculate. See `stoplight_properties_help.txt` for details
`--drop_invalid` by default, invalid smiles will be saved in the output, but listed as invalid and given NA values. Using this drops these (silently) from the output file

# Webserver interface
This repository also contains the code to run a local webserver (or host your own). You can start the server by running `qunicorn wsqi:app` (or using the development flask server by setting the `FLASK_APP` variable: `$env:FLASK_APP = "app"` on windows or `export FLASK_ENV=app` on unix). From that access 127.0.0.1:5000 to view the local server

Thanks to JSME for a free and easy to use molecule editor for webpages
Bienfait, B., Ertl, P. JSME: a free molecule editor in JavaScript. J Cheminform 5, 24 (2013). https://doi.org/10.1186/1758-2946-5-24
