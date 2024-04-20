const moleculeWrapper = document.querySelector('.molecule-wrapper');
const moleculeProperties = document.querySelector('.molecule-properties');
const moleculeSVG = document.getElementById('molecule-svg');
const moleculeSMILES = document.getElementById('molecule_smile_string');
const loadingWrapper = document.querySelector('.loading-wrapper');
const stopLight = document.getElementById("stoplight-svg");
const resDisply = document.getElementById("display-4");

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

export function hideLoadingWrapper() {
    if (!loadingWrapper.className.includes('hidden')) {
        loadingWrapper.classList.add('hidden');
    }
}

export function displayMoleculeCard(moleculeData) {
    hideLoadingWrapper();
    showMoleculeWrapper();
    moleculeProperties.innerHTML = '';
    moleculeSMILES.innerHTML = moleculeData.SMILES;
    moleculeSVG.innerHTML = moleculeData.svg;
    stopLight.innerHTML = moleculeData.stoplight;
    resDisply.innerHTML = `Result: Overall score ${moleculeData.overall}`;

    let groupNum = -1;
    let groupChange = -1;

    for (const [modelName, gnum, classification, confidence, ad, score, color] of moleculeData.molProperties) {

        if (groupNum === -1) {
            groupNum = gnum;
            groupChange = 1;
        }

        if (groupNum !== gnum) {
            groupNum = gnum;
            groupChange = 1;
        }

        if (groupChange === 1) {
            if (groupNum === 0) {
                let propHeader = document.createElement('h5');
                propHeader.classList.add("prop-class-header");
                propHeader.innerHTML = "Molecular Properties";
                moleculeProperties.append(propHeader);
                groupChange = 0;
            }
            if (groupNum === 1) {
                let propHeader = document.createElement('h5');
                propHeader.classList.add("prop-class-header");
                propHeader.innerHTML = "Assay Liabilities";
                moleculeProperties.append(propHeader);
                let conHeader = document.createElement('span');
                conHeader.classList.add("model-preds");
                conHeader.innerHTML = `<strong>Liability Consensus (average probability of liability occuring): ${moleculeData.assayCon}</strong>`;
                moleculeProperties.append(conHeader);
                groupChange = 0;
            }
            if (groupNum === 2) {
                let propHeader = document.createElement('h5');
                propHeader.classList.add("prop-class-header");
                propHeader.innerHTML = "PK Properties";
                moleculeProperties.append(propHeader);
                groupChange = 0;
            }
        }

        let wrapper = document.createElement('div');
        wrapper.className = 'option-item custom-control custom-checkbox mb-3';

        let information = document.createElement('span');
        information.id = modelName;
        information.classList.add('model-preds');

        if (confidence === "") {
            information.innerHTML = `<span>${modelName}: <span style="color:${color};">${classification}</span> | Score: ${score}</span>`;
        } else {
            if (ad === "") {
                if (score === "") {
                    information.innerHTML = `<span>${modelName}: <span style="color:${color};">${classification}</span> | Confidence: ${confidence}</span>`;
                } else {
                    information.innerHTML = `<span>${modelName}: <span style="color:${color};">${classification}</span> | Confidence: ${confidence} | Score: ${score}</span>`;
                }
            } else {
                if (score === "") {
                    information.innerHTML = `<span>${modelName}: <span style="color:${color};">${classification}</span> | Confidence: ${confidence} | AD: ${ad}</span>`;
                } else {
                    information.innerHTML = `<span>${modelName}: <span style="color:${color};">${classification}</span> | Confidence: ${confidence} | AD: ${ad} | Score: ${score}</span>`;
                }
            }
        }

        // let line_break = document.createElement('hr')
        // line_break.classList.add("style1")
        // information.append(line_break)

        wrapper.append(information);
        moleculeProperties.append(wrapper);
    }
}
