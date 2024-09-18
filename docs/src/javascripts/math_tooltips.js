function addTooltips() {
    let refElements = findRefElements();
    refElements.forEach(refElement => addTooltip(refElement));
};


function addTooltip(refElement) {
    let label = getRefLabel(refElement);
    let equation = findEquationByLabel(label);
    if (equation) {
        let tooltip = createTooltip(equation);
        document.getElementsByTagName("p")[0].parentElement.appendChild(tooltip);
        refElement.addEventListener("mouseenter", (ev) => displayTooltipAndUpdatePosition(tooltip, ev.target));
        refElement.addEventListener("mouseleave", () => hideTooltip(tooltip));
    } else {
        console.error("Equation not found with label: " + label);
    } 
};

function findRefElements() {
    let refElements = document.getElementsByClassName("MathJax_ref");
    return Array.from(refElements).map(getParent).filter(isAnchor).filter(isNotDulpicateRefElement);
};
function getParent(element) {
    return element.parentElement;
};
function isAnchor(element) {
    return element.tagName === "A";
};

function isNotDulpicateRefElement(element) {
    return element.parentElement.parentElement.parentElement.tagName !== "MJX-ASSISTIVE-MML";
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
    return labelElement.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement.parentElement;
};


function getRefLabel(refElement) {
    if (refElement.href.at(-1) === "#") {
        console.log(refElement);
    }
    return refElement.href.split("#mjx-eqn%3A").at(-1).replaceAll("%3A", ":");
};


function createTooltip(equationElement) {
    let tooltip = document.createElement("div");
    tooltip.className = "tooltip-math arithmatex";
    tooltip.innerHTML = equationElement.innerHTML;
    removeLabel(tooltip);
    return tooltip;
}
function removeLabel(mathElement) {
    let labelElement = mathElement.children[0].children[0].children[0].children[1];
    mathElement.children[0].children[0].children[0].removeChild(labelElement);
};


function displayTooltipAndUpdatePosition(tooltip, refElement) {
    tooltip.style.display = "block";
    setTooltipPosition(tooltip, refElement);
}
function hideTooltip(tooltip) {
    tooltip.style.display = "none";
}
function setTooltipPosition(tooltip, refElement) {
    let refRect = refElement.getBoundingClientRect();
    let tooltipRect = tooltip.getBoundingClientRect();
    tooltip.style.left = calcLeft(refRect, tooltipRect) + "px";
    tooltip.style.top = calcTop(refRect, tooltipRect) + "px";
}
function calcLeft(refRect, tooltipRect) {
    if (refRect.left + tooltipRect.width > window.innerWidth) {
        return window.innerWidth - tooltipRect.width;
    } else {
        return refRect.left;
    }
}
function calcTop(refRect, tooltipRect) {
    return window.scrollY + refRect.bottom + 10
}


document.addEventListener("MathJaxFinished", addTooltips);