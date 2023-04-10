import { displayError } from './error.js';

const optionsForm = document.getElementById('search-options');
const openButton = document.getElementById('options-dropdown-open');
const closeButton = document.getElementById('options-dropdown-close');

function showOptions() {
    optionsForm.classList.add('search-options-open');

    if (optionsForm.className.includes('search-options-closed')) {
        optionsForm.classList.remove('search-options-closed');
    }

    openButton.classList.add('hidden');
    closeButton.classList.remove('hidden');
}

function hideOptions() {
    optionsForm.classList.add('search-options-closed');

    if (optionsForm.className.includes('search-options-open')) {
        optionsForm.classList.remove('search-options-open');
    }

    closeButton.classList.add('hidden');
    openButton.classList.remove('hidden');
}

openButton.addEventListener('click', showOptions);
closeButton.addEventListener('click', hideOptions);

function addMoleculePropertyOption(name, checked) {
    let wrapper = document.createElement('div');
    wrapper.className = 'option-item custom-control custom-checkbox mb-3';

    let checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.checked = checked === 1;
    checkbox.name = name;
    checkbox.id = name;
    checkbox.classList.add('custom-control-input');

    let label = document.createElement('label');
    label.innerText = name;
    label.htmlFor = name;
    label.classList.add('custom-control-label');

    wrapper.append(checkbox);
    wrapper.append(label);
    optionsForm.append(wrapper);
}

export function displayMoleculePropertyOptions(propNames) {
    let groupOne = document.createElement('p');
    groupOne.innerHTML = "Molecular Properties"
    optionsForm.append(groupOne);
    for (const prop of propNames.mp) {
        addMoleculePropertyOption(prop, 1);
    }

    let groupTwo = document.createElement('p');
    groupTwo.innerHTML = "Assay Liabilities"
    optionsForm.append(groupTwo);
    for (const prop of propNames.al) {
        addMoleculePropertyOption(prop, 0);
    }

    let groupThree = document.createElement('p');
    groupThree.innerHTML = "PK Properties"
    optionsForm.append(groupThree);
    for (const prop of propNames.pk) {
        addMoleculePropertyOption(prop, 0);
    }
}

export function getOptions() {
    const optionsFormData = new FormData(optionsForm);
    let options = {};

    for (const [key, value] of optionsFormData.entries()) {
        // Ticked checkboxes have a value of 'on', need to change it to true
        options[key] = value == 'on' ? true : value;
    }

    return options;
}

fetch('/endpoints', {
    method: 'GET',
})
    .then((response) =>
        response.ok ? response.json() : Promise.reject(response)
    )
    .then((data) => displayMoleculePropertyOptions(data))
    .catch((err) => displayError(err));
