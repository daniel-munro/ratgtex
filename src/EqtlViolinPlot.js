/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */

"use strict";
import { json } from "d3-fetch";
import { createSvg } from "./modules/utils";
import { getGtexUrls, parseDynEqtl } from "./modules/DataParser";
import GroupedViolin from "./modules/GroupedViolin";

export function render(
  par,
  geneId,
  variantId,
  tissueId,
  groupName = undefined,
  urls = getGtexUrls()
) {
  json(
    urls["dyneqtl"] +
      `variantId=${variantId}&geneId=${geneId}&tissueSiteDetailId=${tissueId}`,
    { credentials: "include" }
  ).then(function (jsonResp) {
    const data = parseDynEqtl(jsonResp);

    // construct the dynEqtl data for the three genotypes: ref, het, alt
    par.data = [
      {
        group: groupName || data.tissueSiteDetailId,
        label: data.ref.length > 2 ? "ref" : data.ref,
        size: data.homoRefExp.length,
        values: data.homoRefExp,
      },
      {
        group: groupName || data.tissueSiteDetailId,
        label: data.het.length > 2 ? "het" : data.het,
        size: data.heteroExp.length,
        values: data.heteroExp,
      },
      {
        group: groupName || data.tissueSiteDetailId,
        label: data.alt.length > 2 ? "alt" : data.alt,
        size: data.homoAltExp.length,
        values: data.homoAltExp,
      },
    ];

    // Rendering directly with GroupedViolin (decoupled from GTExViz)
    const margin = {
      top: par.marginTop,
      right: par.marginRight,
      bottom: par.marginBottom,
      left: par.marginLeft,
    };
    const inWidth = par.width - (par.marginLeft + par.marginRight);
    const inHeight = par.height - (par.marginTop + par.marginBottom);

    // Create SVG root group inside container
    const svg = createSvg(par.id, par.width, par.height, margin);

    // Instantiate and render
    const gViolin = new GroupedViolin(par.data);
    gViolin.render(
      svg,
      inWidth,
      inHeight,
      par.xPadding,
      undefined, // xDomain
      [], // yDomain (auto)
      par.yLabel,
      par.showX,
      par.xAngle,
      par.showSubX,
      par.subXAngle,
      par.showWhisker,
      par.showDivider,
      par.showLegend,
      par.showSampleSize,
      par.sortSubX,
      par.showOutliers,
      10 // numPoints: show jittered points if fewer than this
    );

    // Hide the size axis used for sample size labels
    svg
      .selectAll(".violin-size-axis")
      .classed("violin-size-axis-hide", true)
      .classed("violin-size-axis", false);

    // Tooltip for hover details
    const tooltipId = `${par.id}Tooltip`;
    gViolin.createTooltip(tooltipId);
  });
}
