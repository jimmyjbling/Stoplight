import { getOptions } from '/static/options.js';
import { hideMoleculeWrapper, displayMoleculeCard } from "./molecule-card.js";
import { clearErrorMessage, displayError } from "./error.js";

const form = document.getElementById('single-smiles-form');
const smilesInput = document.getElementById('smiles-input');

function resetElements() {
    hideMoleculeWrapper();
    clearErrorMessage();
}

const setVisible = (elementOrSelector, visible) =>
  (typeof elementOrSelector === 'string'
    ? document.querySelector(elementOrSelector)
    : elementOrSelector
  ).style.display = visible ? 'block' : 'none';

form.onsubmit = (event) => {
    event.preventDefault();

    resetElements();

    setVisible('#loading', true);

    fetch('/smiles', {
        method: 'POST',
        body: JSON.stringify({'smiles': smilesInput.value, 'options': getOptions()}),
        headers: {
            'Content-type': 'application/json; charset=UTF-8',
        }

    })
    .then((response) => (response.ok ? response.json() : Promise.reject(response)))
    .then((data) => displayMoleculeCard(data))
    .catch((err) => (displayError(err)))
}