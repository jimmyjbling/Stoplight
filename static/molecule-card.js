const moleculeWrapper = document.querySelector('.molecule-wrapper');
const stoplightWrapper = document.querySelector('.stoplight-wrapper');
const moleculeProperties = document.querySelector('.molecule-properties');
const moleculeSVG = document.getElementById('molecule-svg');
const moleculeStoplight = document.getElementById('stoplight-svg');
const moleculeSMILES = document.getElementById('molecule_smile_string');


function addMolecularPropertyElement(name, value, score, color) {
    let element = document.createElement('span');
    element.class = 'molecular-property';

    element.innerHTML = `${name} <p style="color:${color};">${value} | Score: ${score}</p>`;
    moleculeProperties.append(element);
}

export function showMoleculeWrapper() {
    if (moleculeWrapper.className.includes('hidden')) {
        moleculeWrapper.classList.remove('hidden');
    }
}
export function hideMoleculeWrapper() {
    if (!moleculeWrapper.className.includes('hidden')) {
        moleculeWrapper.classList.add('hidden');
    }
}

export function showStoplightWrapper() {
    if (stoplightWrapper.className.includes('hidden')) {
        stoplightWrapper.classList.remove('hidden');
    }
}
export function hideStoplightWrapper() {
    if (!stoplightWrapper.className.includes('hidden')) {
        stoplightWrapper.classList.add('hidden');
    }
}

export function displayMoleculeCard(moleculeData) {
    showMoleculeWrapper();
    showStoplightWrapper();
    moleculeProperties.innerHTML = '';
    moleculeSMILES.innerHTML = moleculeData.SMILES;
    moleculeSVG.innerHTML = moleculeData.svg;
    moleculeStoplight.innerHTML = moleculeData.stoplight

    for (const [propName, value, score, color] of moleculeData.molProperties) {
        addMolecularPropertyElement(propName, value, score, color);
    }
}