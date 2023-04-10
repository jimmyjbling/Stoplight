import { getOptions } from './options.js';
import { hideMoleculeWrapper, displayMoleculeCard } from "./molecule-card.js";
import { clearErrorMessage, displayError } from "./error.js";

const form = document.getElementById('single-smiles-form');
const smilesInput = document.getElementById('smiles-input');
const loadingWrapper = document.querySelector('.loading-wrapper');

const errorMessage = document.getElementById('error-message-input');
const errorWrapper = document.querySelector('.error-wrapper-input');

function resetElements() {
    hideMoleculeWrapper();
    clearErrorMessage();
}

export function showLoadingWrapper() {
    if (loadingWrapper.className.includes('hidden')) {
        loadingWrapper.classList.remove('hidden');
    }
}

function checkSmiles(data) {
    console.log(data)
    let is_smiles = data.is_smile
    if (is_smiles === 0) {
        errorWrapper.classList.remove('hidden');
        errorMessage.innerHTML = 'Invalid SMILES Input!';
        smilesInput.classList.add('is-invalid');
        form.classList.add('has-danger');
    } else {
        clearErrorMessage();
    }
}

export function checkSmiles2() {
    console.log("check")
}

smilesInput.addEventListener('change', () => {
    clearErrorMessage();
    const res = fetch('/check_smiles', {
        method: 'POST',
        body: JSON.stringify({'smiles': smilesInput.value}),
        headers: {
            'Content-type': 'application/json; charset=UTF-8',
        }
    })
    .then((response) => (response.ok ? response.json() : Promise.reject(response)))
    .then((data) => checkSmiles(data))
    .catch((err) => (displayError(err)))
});

form.onsubmit = (event) => {
    event.preventDefault();

    resetElements();

    showLoadingWrapper();
    loadingWrapper.scrollIntoView({behavior:"smooth"});

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