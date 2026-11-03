/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
/**
 * Create a toolbar
 * This class uses a lot of jQuery for dom element manipulation
 */

import { select } from "d3-selection";
import { downloadSvg } from "./utils.js";

export default class Toolbar {
  constructor(domId, tooltip = undefined, vertical = false) {
    $(`#${domId}`).show(); // if hidden

    // add a new bargroup div to domID with bootstrap button classes
    const btnClasses = vertical
      ? "btn-group-vertical btn-group-sm"
      : "btn-group btn-group-sm";
    this.bar = $("<div/>").addClass(btnClasses).appendTo(`#${domId}`);
    this.buttons = {};
    this.tooltip = tooltip;
  }

  /**
   * Create a download button for SVG
   * @param id {String} the button dom ID
   * @param svgId {String} the SVG dom ID to grab and download
   * @param outfileName {String} the download file name
   * @param cloneId {String} the cloned SVG dom ID
   * @param icon {String} a fontawesome's icon class name
   */
  createDownloadSvgButton(
    id,
    svgId,
    outfileName,
    cloneId,
    icon = "fa-download"
  ) {
    this.createButton(id, icon);
    select(`#${id}`)
      .on("click", () => {
        this.downloadSvg(svgId, outfileName, cloneId);
      })
      .on("mouseover", (event) => {
        this.tooltip?.show("Download", event);
      })
      .on("mouseout", () => {
        this.tooltip?.hide();
      });
  }

  createResetButton(id, callback, icon = "fa-expand-arrows-alt") {
    this.createButton(id, icon);
    select(`#${id}`)
      .on("click", callback)
      .on("mouseover", (event) => {
        this.tooltip?.show("Reset the scales", event);
      })
      .on("mouseout", () => {
        this.tooltip?.hide();
      });
  }

  /**
   * create a button to the toolbar
   * @param id {String} the button's id
   * @param icon {String} a fontawesome icon class
   * Dependencies: Bootstrap, jQuery, Fontawesome
   */
  createButton(id, icon = "fa-download") {
    const $button = $("<a/>")
      .attr("id", id)
      .addClass("btn btn-outline-secondary")
      .appendTo(this.bar);
    $("<i/>").addClass(`fa ${icon}`).appendTo($button);
    this.buttons[id] = $button;
    return $button;
  }

  /**
   * attach a tooltip dom with the toolbar
   * @param tooltip {Tooltip}
   */
  attachTooltip(tooltip) {
    this.tooltip = tooltip;
  }

  /**
   * Download SVG obj
   * @param svgId {String} the SVG dom ID
   * @param fileName {String} the output file name
   * @param cloneId {String} the temporary dom ID to copy the SVG to
   * Dependencies: FileSaver
   */
  downloadSvg(svgId, fileName, cloneId) {
    downloadSvg($(`#${svgId}`), fileName, cloneId);
  }
}
