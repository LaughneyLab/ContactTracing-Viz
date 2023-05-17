// We can't import d3 directly, so we have to use the global variable

window.past_circos_transform = getCircosTransform();

// Use mutation observer
window.circos_observer = new MutationObserver(function (mutations) {
    if (mutations.length > 0) {
        mutation = mutations[mutations.length - 1];
        if (mutation.type === "attributes") {
            past_circos_transform = getCircosTransform();
        }
    }
});

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


function setCircosTransform(transform) {
    // See for more info: https://github.com/nicgirault/circosJS/commit/2f651f5c92ae737c2d036c3b0f14f9441e39fc25
    const circos_top_level = getCircosSvgElement();
    if (circos_top_level !== null) {
        updated = d3.select(circos_top_level)
            .attr("transform", "translate(" + transform.translate + ") scale(" + transform.scale + ")")

        zoom = updated.node().parentNode.__zoom;
        if (zoom !== undefined) {
            zoom.x = transform.translate[0];
            zoom.y = transform.translate[1];
            zoom.k = transform.scale;
        }

        past_circos_transform = transform;
    }
}

function resetCircosTransform() {
    const circos_top_level = getCircosSvgElement().parentElement;
    setCircosTransform({ translate: [circos_top_level.getAttribute("width")/2, circos_top_level.getAttribute("height")/2], scale: 1.0 });
}

// Hook the circos plot
function moveCircosTooltip() {
    let tooltip = window.document.querySelector(".circos-tooltip");
    const circos_top_level = getCircosSvgElement();
    if (tooltip !== null && circos_top_level !== null) {
        // Move the tooltip to be a child of the circos div
        circos_container = circos_top_level.parentElement.parentElement;
        tooltip.parentElement.removeChild(tooltip);
        circos_container.appendChild(tooltip);
    }
}

function checkPartners(partnersObj, celltype, target) {
    if (partnersObj === undefined) {
        return false;
    }
    let partners = partnersObj[celltype];
    if (partners === undefined) {
        return false;
    }
    return partners.includes(target);
}

function injectHoverEffects() {
    // Highlight relevant connections on hover
    let svg = d3.select("#svg-child");
    // Add effects by applying the not-hovered class to all non-hovered elements
    svg.selectAll(".chord").on("mouseover", function (d) {  // Connections
        svg.classed("moused-over", true);
        // Get the properties of the hovered chord
        let source = d.source_celltype;
        let target = d.target_celltype;
        let ligand = d.ligand;
        let receptor = d.receptor;

        svg.selectAll("path").filter(function (d) {  // Select tracks that don't contain the hovered chord properties
            // Ignore elements with undefined properties
            if (d.celltype === undefined || d.target === undefined) {
                return false;
            }
            return !((d.celltype === source && d.target === ligand) || (d.celltype === target && d.target === receptor));
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".cs-layout > *").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return d.celltype !== source && d.celltype !== target;
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".chord").filter(function (d) {  // Select chords that don't contain the hovered chord properties
            return !(d.source_celltype === source && d.ligand === ligand && d.target_celltype === target && d.receptor === receptor);
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".block > g > text").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return !((d.celltype === source && d.value === ligand) || (d.celltype === target && d.value === receptor));
        }).classed("not-hovered", true);  // Apply the not-hovered class
    }).on("mouseout", function (d) {
        svg.selectAll(".not-hovered").classed("not-hovered", false);  // Remove the not-hovered class
        svg.classed("moused-over", false);
    });

    svg.selectAll("path").filter(function(d) {
        return d.celltype !== undefined && d.target !== undefined;
    }).on("mouseover", function (d) {  // Tracks
        svg.classed("moused-over", true);
        // Get the properties of the hovered chord
        let celltype = d.celltype;
        let target = d.target;
        let partners = d.partners;

        svg.selectAll("path").filter(function (d) {  // Select tracks that don't contain the hovered chord properties
            // Ignore elements with undefined properties
            if (d.celltype === undefined || d.target === undefined) {
                return false;
            }
            return !((d.celltype === celltype && d.target === target) || checkPartners(partners, d.celltype, d.target));
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".cs-layout > *").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return d.celltype !== celltype && partners[d.celltype] === undefined;
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".chord").filter(function (d) {  // Select chords that don't contain the hovered chord properties
            return !((d.source_celltype === celltype && d.ligand === target) || (d.target_celltype === celltype && d.receptor === target));
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".block > g > text").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return !((d.celltype === celltype && d.value === target) || checkPartners(partners, d.celltype, d.value));
        }).classed("not-hovered", true);  // Apply the not-hovered class
    }).on("mouseout", function (d) {
        svg.selectAll(".not-hovered").classed("not-hovered", false);  // Remove the not-hovered class
        svg.classed("moused-over", false);
    });

    svg.selectAll(".cs-layout > *").on("mouseover", function (d) {  // Cell type labels
        svg.classed("moused-over", true);
        // Get the properties of the hovered chord
        let celltype = d.celltype;
        let receiving_celltypes = d.receiving_celltypes;
        let sending_celltypes = d.sending_celltypes;

        svg.selectAll("path").filter(function (d) {  // Select tracks that don't contain the hovered chord properties
            // Ignore elements with undefined properties
            if (d.celltype === undefined || d.target === undefined) {
                return false;
            }
            return d.celltype !== celltype && !receiving_celltypes.includes(d.celltype) && !sending_celltypes.includes(d.celltype);
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".cs-layout > *").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return d.celltype !== celltype && !receiving_celltypes.includes(d.celltype) && !sending_celltypes.includes(d.celltype);
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".chord").filter(function (d) {  // Select chords that don't contain the hovered chord properties
            return !((d.source_celltype === celltype) || (d.target_celltype === celltype));
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".block > g > text").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return d.celltype !== celltype && d.partners[celltype] === undefined && d.partners[celltype] === undefined;
        }).classed("not-hovered", true);  // Apply the not-hovered class
    }).on("mouseout", function (d) {
        svg.selectAll(".not-hovered").classed("not-hovered", false);  // Remove the not-hovered class
        svg.classed("moused-over", false);
    });

    svg.selectAll(".block > g > text").on("mouseover", function (d) {  // Cell type labels
        svg.classed("moused-over", true);
        // Get the properties of the hovered chord
        let celltype = d.celltype;
        let target = d.value;
        let partners = d.partners;

        svg.selectAll("path").filter(function (d) {  // Select tracks that don't contain the hovered chord properties
            // Ignore elements with undefined properties
            if (d.celltype === undefined || d.target === undefined) {
                return false;
            }
            return !((d.celltype === celltype && d.target === target) || checkPartners(partners, d.celltype, d.target));
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".cs-layout > *").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return d.celltype !== celltype && partners[d.celltype] === undefined;
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".chord").filter(function (d) {  // Select chords that don't contain the hovered chord properties
            return !((d.source_celltype === celltype && d.ligand === target) || (d.target_celltype === celltype && d.receptor === target));
        }).classed("not-hovered", true);  // Apply the not-hovered class

        svg.selectAll(".block > g > text").filter(function (d) {  // Select cell type labels that don't contain the hovered chord properties
            return !((d.celltype === celltype && d.value === target) || checkPartners(partners, d.celltype, d.value));
        }).classed("not-hovered", true);  // Apply the not-hovered class
    }).on("mouseout", function (d) {
        svg.selectAll(".not-hovered").classed("not-hovered", false);  // Remove the not-hovered class
        svg.classed("moused-over", false);
    });

    // On double click, undo highlighting
    svg.on("dblclick", function (d) {
        svg.selectAll('.not-hovered').classed("not-hovered", false)
        svg.classed("moused-over", false);
    });
}

function circosInjection() {
    moveCircosTooltip();

    injectHoverEffects();

    circos_observer.disconnect();

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

function zoomInCircos() {
    zoomCircos(1.1);
}

function zoomOutCircos() {
    zoomCircos(1 / 1.1);
}

// Expose the functions globally
window.circosInjection = circosInjection;
window.downloadCircosSvg = downloadCircosSvg;
window.resetCircosTransform = resetCircosTransform;
window.zoomInCircos = zoomInCircos;
window.zoomOutCircos = zoomOutCircos;
window.moveCircosTooltip = moveCircosTooltip;
