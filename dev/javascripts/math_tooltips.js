function addTooltips() {
    let refElements = findRefElements();
    refElements.forEach(refElement => addTooltip(refElement));
};


function addTooltip(refElement) {
    let label = getRefLabel(refElement);
    let equation = findEquationByLabel(label);
    if (equation) {
        let tooltip = createTooltip(equation);
        refElement.appendChild(tooltip);
        refElement.addEventListener("mouseover", () => displayTooltip(tooltip));
        refElement.addEventListener("mouseout", () => hideTooltip(tooltip));
    } else {
        console.error("Equation not found with label: " + label);
    } 
};

function findRefElements() {
    let refElements = document.getElementsByClassName("MathJax_ref");
    return Array.from(refElements).map(getParent).filter(isAnchor);
};
function getParent(element) {
    return element.parentElement;
};
function isAnchor(element) {
    return element.tagName === "A";
};


function findEquationByLabel(label) {
    let labelElement = getLabelElement(label);
    if (labelElement) {
        return getMathElement(labelElement);
    } else {
        console.error("Equation with label " + label + " not found");
    }
};
function getLabelElement(label) {
    let labelId = "mjx-eqn:" + label;
    return document.getElementById(labelId);
};
function getMathElement(labelElement) {
    return labelElement.parentElement.parentElement.parentElement.parentElement.children[0].children[0];
};


function getRefLabel(refElement) {
    if (refElement.href.at(-1) === "#") {
        console.log(refElement);
    }
    return refElement.href.split("#mjx-eqn%3A").at(-1).replaceAll("%3A", ":");
};


function createTooltip(equation) {
    let tooltip = document.createElement("div");
    tooltip.className = "tooltip-math";
    tooltip.innerHTML = equation.outerHTML;
    return tooltip;
}


function displayTooltip(tooltip) {
    tooltip.style.display = "block";
}
function hideTooltip(tooltip) {
    tooltip.style.display = "none";
}


document.addEventListener("MathJaxFinished", addTooltips);