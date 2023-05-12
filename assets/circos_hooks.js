// We can't import d3 directly, so we have to use the global variable
d3 = window.d3;

// Retrieve the circos transform from the DOM
function getCircosSvgElement() {
    let doc = window.document;
    return doc.querySelector("#Circos-container > #svg-child > g");
}

function getCircosTransform() {
    let scale = 1.0;
    let translate = [0, 0];
    const circos_top_level = getCircosSvgElement();
    if (circos_top_level !== null) {
        let circos_transform = circos_top_level.getAttribute("transform")
        // Try to parse the transform and scale
        if (circos_transform !== null) {
            let split = circos_transform.split("scale(");
            translate = split[0].replace("translate(", "").replace(")", "").split(",").map(function (x) {
                return parseFloat(x);
            });

            if (split.length > 1) {
                scale = parseFloat(split[1].replace(")", ""));
            }
        }
    }
    return { translate: translate, scale: scale };
}

past_circos_transform = getCircosTransform();

// Use mutation observer
circos_observer = new MutationObserver(function (mutations) {
    if (mutations.length > 0) {
        if (mutation.type === "attributes") {
            mutation = mutations[mutations.length - 1];
            past_circos_transform = getCircosTransform();
        }
    }
});

function setCircosTransform(transform) {
    // See for more info: https://github.com/nicgirault/circosJS/commit/2f651f5c92ae737c2d036c3b0f14f9441e39fc25
    const circos_top_level = getCircosSvgElement();
    if (circos_top_level !== null) {
        // Use d3 to zoom and transform
        d3.select(circos_top_level)
            .attr("transform", "translate(" + transform.translate + ") scale(" + transform.scale + ")");
        past_circos_transform = transform;
    }
}

function resetCircosTransform() {
    const circos_top_level = getCircosSvgElement().parentElement;
    setCircosTransform({ translate: [circos_top_level.getAttribute("width")/2, circos_top_level.getAttribute("height")/2], scale: 1.0 });
}

// Hook the circos plot
function circosInjection() {
    circos_observer.disconnect();

    // Global variable to store previous transformation
    setCircosTransform(past_circos_transform);

    // Add listener for when the transform attribute changes
    const circos_top_level = getCircosSvgElement();
    if (circos_top_level !== null) {
        circos_observer.observe(circos_top_level, { attributes: true });
    }
}

function downloadCircosSvg() {
    let svg = getCircosSvgElement().parentElement;
    svgExport.downloadSvg(
        svg,
        "circos",
        {
            useCSS: false
        }
    );
}

function zoomCircos(zoom) {
    let transform = getCircosTransform();
    transform.scale *= zoom;
    setCircosTransform(transform);
}

CIRCOS_ZOOM_FACTOR = 1.1;

function zoomInCircos() {
    zoomCircos(CIRCOS_ZOOM_FACTOR);
}

function zoomOutCircos() {
    zoomCircos(1 / CIRCOS_ZOOM_FACTOR);
}
