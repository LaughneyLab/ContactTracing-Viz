function getSvgElements() {
    let doc = window.document;
    return doc.querySelectorAll("div.js-plotly-plot > div > div > *");
}

function plotlyInjection(heightRatio) {
    // Add listener for when the transform attribute changes
    return;  // Current unneeded with current CSS settings
    const baseSize = 75;
    const width = baseSize;
    const height = baseSize * heightRatio;
    const viewBoxWidth = 1275;
    const viewBoxHeight = 1275 * heightRatio;
    const svgs = getSvgElements();
    if (svgs.length === 0) {
        // Plotly hasn't rendered the SVG yet, try again in 10ms
        setTimeout(plotlyInjection, 10, heightRatio);
        return;
    }
    svgs.forEach(function (svg) {
        // If icon, adjust viewbox to translate to center (this is for the modebar)
        svg.setAttribute('viewBox', '0 0 ' + viewBoxWidth + ' ' + viewBoxHeight);
        svg.setAttribute('width', width + 'vw');
        svg.setAttribute('height', height + 'vw');
        svg.setAttribute('transform', 'translate(0, 0) scale(1)');
    });
}

// Expose the functions globally
window.plotlyInjection = plotlyInjection;
