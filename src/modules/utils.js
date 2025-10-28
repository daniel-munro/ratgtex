/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
/**
 * Creates an SVG
 * @param id {String} a DOM element ID that starts with a "#"
 * @param width {Numeric}
 * @param height {Numeric}
 * @param margin {Object} with two attributes: width and height
 * @return {Selection} the d3 selection object of the SVG
 */

import { select } from "d3-selection";

export function checkDomId(id) {
  // test input params
  if ($(`#${id}`).length == 0) {
    let error = `Input Error: DOM ID ${id} is not found.`;
    //alert(error);
    console.warn(error);
    throw error;
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
  svgId = undefined,
  position = "relative"
) {
  checkDomId(id);
  if (svgId === undefined) svgId = `${id}-svg`;
  return (
    select("#" + id)
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("id", svgId)
      // .style("position", position)
      .append("g")
      .attr("transform", `translate(${margin.left}, ${margin.top})`)
  );
}

/**
 *
 * @param svgObj
 * @param downloadFileName {String}
 * @param tempDownloadDivId {String}
 */
export function downloadSvg(svgObj, downloadFileName, tempDownloadDivId) {
  console.log(svgObj);
  var $svgCopy = svgObj
    .clone()
    .attr("version", "1.1")
    .attr("xmlns", "http://www.w3.org/2000/svg");

  // parse and add the CSS styling used by the SVG
  var styles = parseCssStyles(svgObj.get());
  $svgCopy.prepend(styles);

  $("#" + tempDownloadDivId)
    .html("")
    .hide();
  var svgHtml = $("#" + tempDownloadDivId)
    .append($svgCopy)
    .html();

  var svgBlob = new Blob([svgHtml], { type: "image/svg+xml" });
  saveAs(svgBlob, downloadFileName);

  // clear the temp download div
  $("#" + tempDownloadDivId)
    .html("")
    .hide();
}
/**
 * A function for parsing the CSS style sheet and including the style properties in the downloadable SVG.
 * @param dom
 * @returns {Element}
 */
export function parseCssStyles(dom) {
  var used = "";
  var sheets = document.styleSheets;

  for (var i = 0; i < sheets.length; i++) {
    // TODO: walk through this block of code

    try {
      if (sheets[i].cssRules == null) continue;
      var rules = sheets[i].cssRules;

      for (var j = 0; j < rules.length; j++) {
        var rule = rules[j];
        if (typeof rule.style != "undefined") {
          var elems;
          // removing any selector text including svg element ID -- dom already selects for that
          var selector =
            rule.selectorText === undefined
              ? rule.selectorText
              : rule.selectorText.replace(`#${dom[0].id} `, "");
          //Some selectors won't work, and most of these don't matter.
          try {
            elems = $(dom).find(selector);
          } catch (e) {
            elems = [];
          }

          if (elems.length > 0) {
            used += rule.selectorText + " { " + rule.style.cssText + " }\n";
          }
        }
      }
    } catch (e) {
      // In Firefox, if stylesheet originates from a diff domain,
      // trying to access the cssRules will throw a SecurityError.
      // Hence, we must use a try/catch to handle this in Firefox
      if (e.name !== "SecurityError") throw e;
      continue;
    }
  }

  var s = document.createElement("style");
  s.setAttribute("type", "text/css");
  s.innerHTML = "<![CDATA[\n" + used + "\n]]>";

  return s;
}
