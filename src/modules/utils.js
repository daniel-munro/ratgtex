/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
import { select } from "d3-selection";

export function checkDomId(id) {
  // test input params
  if ($(`#${id}`).length == 0) {
    throw new Error(`Input Error: DOM ID ${id} is not found.`);
  }
}

/**
 * Create an SVG D3 object
 * @param id {String} the parent dom ID
 * @param width {Numeric}: the outer width
 * @param height {Numeric}: the outer height
 * @param margin {Object} with attr: left, top
 * @param svgId {String}
 * @returns {*}
 */
export function createSvg(
  id,
  width,
  height,
  margin,
  svgId = undefined
) {
  checkDomId(id);
  if (svgId === undefined) svgId = `${id}-svg`;
  return select("#" + id)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("id", svgId)
    .append("g")
    .attr("transform", `translate(${margin.left}, ${margin.top})`);
}

/**
 *
 * @param svgObj
 * @param downloadFileName {String}
 */
export function downloadSvg(svgObj, downloadFileName) {
  const original = svgObj.get(0);
  if (!original) throw new Error("SVG export: source element not found.");
  const copy = original.cloneNode(true);
  copy.setAttribute("version", "1.1");
  copy.setAttributeNS("http://www.w3.org/2000/xmlns/", "xmlns", "http://www.w3.org/2000/svg");
  copy.prepend(parseCssStyles(original));

  const svg = new XMLSerializer().serializeToString(copy);
  saveAs(new Blob([svg], { type: "image/svg+xml" }), downloadFileName);
}

/** Collect accessible stylesheet rules that apply to the SVG or its children. */
export function parseCssStyles(dom) {
  const svg = dom.nodeType ? dom : dom[0];
  const collect = (rules) => {
    let css = "";
    for (const rule of rules) {
      if (rule.selectorText && rule.style) {
        try {
          if (svg.matches(rule.selectorText) || svg.querySelector(rule.selectorText)) {
            css += rule.cssText + "\n";
          }
        } catch {
          // Ignore selectors that the browser cannot query (e.g. pseudo-elements).
        }
      } else if (rule.cssRules) {
        const nested = collect(rule.cssRules);
        if (nested) css += rule.cssText.slice(0, rule.cssText.indexOf("{")) + "{\n" + nested + "}\n";
      }
    }
    return css;
  };

  let css = "";
  for (const sheet of document.styleSheets) {
    try {
      if (sheet.cssRules) css += collect(sheet.cssRules);
    } catch (error) {
      // Browsers deny access to rules in cross-origin stylesheets.
      if (error.name !== "SecurityError") throw error;
    }
  }
  const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
  style.setAttribute("type", "text/css");
  style.textContent = css;
  return style;
}
