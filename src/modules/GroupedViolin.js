/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
/*
Input data structure: a list of data object with the following structure:
[
    {
        group: "group1"
        label: "dataset 1",
        values: [a list of numerical values]
     },
     {
        group: "group1"
        label: "dataset 2",
        values: [a list of numerical values]
     },
     {
        group: "group2"
        label: "dataset 3",
        values: [a list of numerical values]
     }
]
*/

import { extent, median, ascending, quantile, max, min } from "d3-array";
import { nest } from "d3-collection";
import { scaleBand, scaleLinear } from "d3-scale";
import { area } from "d3-shape";
import { axisBottom, axisLeft } from "d3-axis";
import { select } from "d3-selection";
import { randomNormal } from "d3-random";

import { kernelDensityEstimator, kernel, kernelBandwidth } from "./kde.js";
import Tooltip from "./Tooltip.js";
import Toolbar from "./Toolbar.js";

export default class GroupedViolin {
  /**
   * constructor for GroupedViolin
   * @param data {List}: a list of objects with attributes: group: {String}, label: {String}, values: {List} of numerical values, size: integer, optional
   * @param groupInfo {Dictionary}: metadata of the group, indexed by group ID
   */
  constructor(data, groupInfo = {}) {
    this._sanityCheck(data);
    this.data = data;
    this.groupInfo = groupInfo;
    this.toolbar = undefined;
    this.tooltip = undefined;
  }

  /**
   *
   * @param dom {DOM} the SVG dom object to append the violin plot to
   * @param width {Float}
   * @param height {Float}
   * @param xPadding {Float} padding of the x axis
   * @param xDomain {List} the order of X groups
   * @param yDomain  {List} the min and max values of the y domain. Set to empty array or undefined to auto-calculate.
   * @param yLabel {String}
   * @param showX
   * @param showSubX
   * @param subXAngle
   * @param showWhisker
   * @param showDivider
   * @param showLegend
   * @param showSize
   * @param sortSubX
   * @param showOutliers
   * @param numPoints {Integer} Number of points required to render violin plot. Displays points if n < numPoints.
   */

  render(
    dom,
    width = 500,
    height = 357,
    xPadding = 0.05,
    xDomain = undefined,
    yDomain = [-3, 3],
    yLabel = "Y axis",
    showX = true,
    xAngle = 30,
    showSubX = true,
    subXAngle = 0,
    showWhisker = false,
    showDivider = false,
    showLegend = false,
    showSize = false,
    sortSubX = false,
    showOutliers = false,
    numPoints = 0
  ) {
    // define the reset for this plot
    this.reset = () => {
      dom.selectAll("*").remove();
      this.render(
        dom,
        width,
        height,
        xPadding,
        xDomain,
        yDomain,
        yLabel,
        showX,
        xAngle,
        showSubX,
        subXAngle,
        showWhisker,
        showDivider,
        showLegend,
        showSize,
        sortSubX,
        showOutliers,
        numPoints
      );
    };

    this.updateYScale = function (yLabel = undefined, yScale = undefined) {
      dom.selectAll("*").remove();
      if (yScale === undefined) yScale = [];
      this.render(
        dom,
        width,
        height,
        xPadding,
        xDomain,
        yScale,
        yLabel,
        showX,
        xAngle,
        showSubX,
        subXAngle,
        showWhisker,
        showDivider,
        showLegend,
        showSize,
        sortSubX,
        showOutliers,
        numPoints
      );
    };

    this.updateXScale = function (xDomain = undefined) {
      if (xDomain === undefined)
        console.error("updateXScale called without new X domain");
      dom.selectAll("*").remove();
      this.render(
        dom,
        width,
        height,
        xPadding,
        xDomain,
        [],
        yLabel,
        showX,
        xAngle,
        showSubX,
        subXAngle,
        showWhisker,
        showDivider,
        showLegend,
        showSize,
        sortSubX,
        showOutliers,
        numPoints
      );
    };

    // defines the X, subX, Y, Z scales
    if (yDomain === undefined || 0 == yDomain.length) {
      let allV = [];
      this.data.forEach((d) => (allV = allV.concat(d.values)));
      yDomain = extent(allV);
    }

    // re-organized this.data indexed by groups
    this.groups = nest()
      .key((d) => d.group)
      .entries(this.data);

    this.scale = {
      x: scaleBand()
        .rangeRound([0, width])
        .domain(xDomain || this.groups.map((d) => d.key))
        .paddingInner(xPadding),
      subx: scaleBand(),
      y: scaleLinear().rangeRound([height, 0]).domain(yDomain),
      z: scaleLinear(), // this is the violin width, the domain and range are determined later individually for each violin
    };

    // for each group
    this.groups.forEach((g) => {
      let group = g.key;
      let entries = g.values;
      let info = this.groupInfo[group]; // optional
      g.index = this.scale.x.domain().indexOf(group);

      if (info !== undefined) {
        // renders group info such as p-value, group name
        const groupInfoDom = dom.append("g");
        const groupLabels = groupInfoDom
          .selectAll(".violin-group-label")
          .data(["pvalue"]);
        groupLabels
          .enter()
          .append("text") // Code review: consider moving this part to the eQTL dashboard
          .attr("x", 0)
          .attr("y", 0)
          .attr("class", "violin-group-label")
          .attr("text-anchor", "center")
          .attr("fill", (d) => {
            // console.log(info['pvalueThreshold']);
            return d == "pvalue" &&
              parseFloat(info[d]) <= parseFloat(info["pvalueThreshold"])
              ? "orangered"
              : "SlateGray";
          })
          .attr("transform", (d, i) => {
            let x = this.scale.x(group) + this.scale.x.bandwidth() / 2;
            let y = this.scale.y(yDomain[0]) + 50; // todo: avoid hard-coded values
            return `translate(${x}, ${y})`;
          })
          // .text((d) => `${d}: ${parseFloat(parseFloat(info[d]).toPrecision(3)).toExponential()}`);
          // .text((d) => `${d}: ${d == 'pvalue' && info[d] === null ? 'NS' : parseFloat(parseFloat(info[d]).toPrecision(3)).toExponential()}`);
          .text((d) =>
            d == "pvalue" && info[d] === null
              ? ""
              : `${d}: ${parseFloat(
                  parseFloat(info[d]).toPrecision(3)
                ).toExponential()}`
          );
      }

      // defines the this.scale.subx based on this.scale.x
      if (sortSubX) {
        entries.sort((a, b) => {
          if (a.label < b.label) return -1;
          else if (a.label > b.label) return 1;
          return 0;
        });
      }

      this.scale.subx
        .domain(entries.map((d) => d.label))
        .rangeRound([
          this.scale.x(group),
          this.scale.x(group) + this.scale.x.bandwidth(),
        ]);

      entries.forEach((entry) => {
        if (0 == entry.values.length) return; // no further rendering if this group has no entries
        entry.values = entry.values.sort(ascending);
        this._drawViolin(
          dom,
          entry,
          showWhisker,
          g.index,
          showOutliers,
          numPoints
        );
      });

      // adds the sub-x axis if there are more than one entries
      var buffer = 18;
      if (showSize) {
        var labels = entries.map((d) => {
          return d.size == undefined ? "(0)" : `(${d.size || 0})`;
        });
        const sizeScale = scaleBand()
          .domain(entries.map((d, i) => i))
          .rangeRound([
            this.scale.x(group),
            this.scale.x(group) + this.scale.x.bandwidth(),
          ]);
        const sizexG = dom
          .append("g")
          .attr("class", "violin-size-axis")
          .attr("transform", `translate(0, ${height + buffer})`)
          .call(axisBottom(sizeScale).tickFormat((i) => labels[i]));
      }

      if (showSubX) {
        var buffer = 5;
        const subxG = dom
          .append("g")
          .attr("class", "violin-sub-axis")
          .attr("transform", `translate(0, ${height + buffer})`)
          .call(axisBottom(this.scale.subx));

        if (subXAngle > 0) {
          subxG
            .selectAll("text")
            .attr("text-anchor", "start")
            .attr("transform", `rotate(${subXAngle}, 2, 10)`);
        }
      }
    });

    // renders the x axis
    let buffer = showSubX ? 55 : 0; // Code review: hard-coded values
    this.xAxis = showX
      ? axisBottom(this.scale.x)
      : axisBottom(this.scale.x).tickFormat("");
    dom
      .append("g")
      .attr("class", "violin-x-axis axis--x")
      .attr("transform", `translate(0, ${height + buffer})`)
      .call(this.xAxis) // set tickFormat("") to show tick marks without text labels
      .selectAll("text")
      .attr("text-anchor", "start")
      .attr("transform", `rotate(${xAngle}, -10, 10)`);

    // adds the y Axis
    buffer = 5;
    this.yAxis = axisLeft(this.scale.y).tickValues(this.scale.y.ticks(5));
    dom
      .append("g")
      .attr("class", "violin-y-axis axis--y")
      .attr("transform", `translate(-${buffer}, 0)`)
      .call(this.yAxis);

    // adds the text label for the y axis
    dom
      .append("text")
      .attr("class", "violin-axis-label")
      .attr("text-anchor", "middle")
      .attr(
        "transform",
        `translate(-${
          buffer * 2 + dom.select(".violin-y-axis").node().getBBox().width
        }, ${
          this.scale.y.range()[0] +
          (this.scale.y.range()[1] - this.scale.y.range()[0]) / 2
        }) rotate(-90)`
      )
      .text(yLabel);

    // plot mouse events
    dom.on("mouseout", () => {
      if (this.tooltip !== undefined) this.tooltip.hide();
    });

    // add group dividers
    if (showDivider) {
      this._addGroupDivider(dom);
    }

    // add color legend
    if (showLegend) {
      const legendG = dom
        .append("g")
        .attr("id", "violinLegend")
        .attr("transform", `translate(0, 0)`);

      legendG
        .append("rect")
        .attr("x", this.scale.x.range()[0])
        .attr("y", -35)
        .attr("width", 60 * this.groups[0].values.length + 10)
        .attr("height", 24)
        .style("fill", "none")
        .style("stroke", "silver");

      const legends = legendG
        .selectAll(".violin-legend")
        .data(this.groups[0].values);

      const g = legends.enter().append("g").classed("violin-legend", true);
      const w = 10;
      g.append("rect")
        .attr("x", (d, i) => 5 + 60 * i + this.scale.x.range()[0])
        .attr("y", -28)
        .attr("width", w)
        .attr("height", w)
        .style("fill", (d) => d.color);

      g.append("text")
        .attr("class", "violin-legend-text")
        .text((d) => d.label)
        .attr("x", (d, i) => 17 + 60 * i + this.scale.x.range()[0])
        .attr("y", -20);
    }
  }

  addPlotTitle(dom, title) {
    let x = (this.scale.x.range()[1] - this.scale.x.range()[0] + 1) / 2;
    let y = this.scale.y.range()[1] - 10;
    dom
      .append("text")
      .attr("class", "violin-title")
      .attr("text-anchor", "middle")
      .attr("transform", `translate(${x}, ${y})`)
      .text(title);
  }

  /**
   * Create the tooltip object
   * @param domId {String} the tooltip's dom ID
   * @returns {Tooltip}
   */
  createTooltip(domId) {
    if ($(`#${domId}`).length == 0)
      $("<div/>").attr("id", domId).appendTo($("body"));
    this.tooltip = new Tooltip(domId);
    select(`#${domId}`).classed("violin-tooltip", true);
    return this.tooltip;
  }

  /**
   * Create the toolbar panel
   * @param domId {String} the toolbar's dom ID
   * @param tooltip {Tooltip}
   * @returns {Toolbar}
   */

  createToolbar(domId, tooltip) {
    // if (tooltip === undefined) tooltip = this.createTooltip(domId);
    this.toolbar = new Toolbar(domId, tooltip);
    return this.toolbar;
  }

  /**
   * render the violin and box plots
   * @param dom {D3 DOM}
   * @param entry {Object} with attrs: values, label
   * @param showWhisker {Boolean}
   * @param gIndex
   * @param showOutliers {Boolean}
   * @param showPoints {Boolean}
   * @private
   */
  _drawViolin(dom, entry, showWhisker, gIndex, showOutliers, numPoints) {
    // generate the vertices for the violin path use a kde
    let kde = kernelDensityEstimator(
      kernel.gaussian,
      this.scale.y.ticks(100), // use up to 100 vertices along the Y axis (to create the violin path)
      kernelBandwidth.nrd(entry.values) // estimate the bandwidth based on the data
    );
    const eDomain = extent(entry.values); // get the max and min in entry.values
    const vertices = kde(entry.values).filter(
      (d) => d[0] >= eDomain[0] && d[0] <= eDomain[1]
    ); // filter the vertices that aren't in the entry.values

    const violinG = dom
      .append("g")
      .attr("id", `violin${gIndex}-${entry.label}`)
      .attr("class", "violin-g")
      .datum(entry);

    // violin plot and box can only be drawn when vertices exist and there are no NaN points
    if (
      entry.values.length >= numPoints &&
      vertices.length &&
      this._validVertices(vertices)
    ) {
      // define the z scale -- the violin width
      let zMax = max(vertices, (d) => Math.abs(d[1])); // find the abs(value) in entry.values
      this.scale.z
        .domain([-zMax, zMax])
        .range([
          this.scale.subx(entry.label),
          this.scale.subx(entry.label) + this.scale.subx.bandwidth(),
        ]);

      // visual rendering
      let violin = area()
        .x0((d) => this.scale.z(d[1]))
        .x1((d) => this.scale.z(-d[1]))
        .y((d) => this.scale.y(d[0]));
      const vPath = violinG
        .append("path")
        .datum(vertices)
        .attr("d", violin)
        .classed("violin", true)
        .style("fill", () => {
          if (entry.color !== undefined) return entry.color;
          // alternate the odd and even colors, maybe we don't want this feature
          if (gIndex % 2 == 0) return "#90c1c1";
          return "#94a8b8";
        });

      // boxplot
      const q1 = quantile(entry.values, 0.25);
      const q3 = quantile(entry.values, 0.75);
      const z = this.scale.z.domain()[1] / 3;

      if (showWhisker) {
        // the upper and lower limits of entry.values
        const iqr = Math.abs(q3 - q1);
        const upper = max(entry.values.filter((d) => d <= q3 + iqr * 1.5));
        const lower = min(entry.values.filter((d) => d >= q1 - iqr * 1.5));
        dom
          .append("line")
          .classed("whisker", true)
          .attr("x1", this.scale.z(0))
          .attr("x2", this.scale.z(0))
          .attr("y1", this.scale.y(upper))
          .attr("y2", this.scale.y(lower))
          .style("stroke", "#fff");
      }

      // interquartile range
      violinG
        .append("rect")
        .attr("x", this.scale.z(-z))
        .attr("y", this.scale.y(q3))
        .attr("width", Math.abs(this.scale.z(-z) - this.scale.z(z)))
        .attr("height", Math.abs(this.scale.y(q3) - this.scale.y(q1)))
        .attr("class", "violin-ir");

      // median
      const med = median(entry.values);
      violinG
        .append("line") // the median line
        .attr("x1", this.scale.z(-z))
        .attr("x2", this.scale.z(z))
        .attr("y1", this.scale.y(med))
        .attr("y2", this.scale.y(med))
        .style("stroke", "#f88")
        .attr("class", "violin-median");

      // outliers
      if (showOutliers) {
        let jitter = randomNormal(0, z / 2);
        const iqr = Math.abs(q3 - q1);
        const upper = max(entry.values.filter((d) => d <= q3 + iqr * 1.5));
        const lower = min(entry.values.filter((d) => d >= q1 - iqr * 1.5));
        const outliers = entry.values.filter((d) => d < lower || d > upper);
        violinG
          .append("g")
          .attr("class", "violin-outliers")
          .selectAll("circle")
          .data(outliers)
          .enter()
          .append("circle")
          // .attr("cx", ()=>this.scale.z(0))
          .attr("cx", () => this.scale.z(jitter()))
          .attr("cy", (d) => this.scale.y(d))
          .attr("r", 2);
      }

      // mouse events
      violinG.on("mouseover", (event) => {
        vPath.classed("highlighted", true);
        if (this.tooltip === undefined)
          console.warn("GroupViolin Warning: tooltip not defined");
        else {
          this.tooltip.show(
            entry.group +
              "<br/>" +
              entry.label +
              "<br/>" +
              "Median: " +
              med.toPrecision(4) +
              "<br/>",
            event
          );
        }
      });
      violinG.on("mouseout", () => {
        vPath.classed("highlighted", false);
      });
    } else if (numPoints) {
      // define the z scale -- the violin width
      let zMax = max(entry.values, (d) => Math.abs(d)); // find the abs(value) in entry.values
      this.scale.z
        .domain([-zMax, zMax])
        .range([
          this.scale.subx(entry.label),
          this.scale.subx(entry.label) + this.scale.subx.bandwidth(),
        ]);
      const z = this.scale.z.domain()[1] / 3;
      let jitter = randomNormal(0, z / 2);

      violinG
        .append("g")
        .attr("class", "violin-points")
        .selectAll("circle")
        .data(entry.values)
        .enter()
        .append("circle")
        .attr("cx", () => this.scale.z(jitter()))
        .attr("cy", (d) => this.scale.y(d))
        .attr("r", 2);
    }
  }

  _sanityCheck(data) {
    const attr = ["group", "label", "values"];

    data.forEach((d) => {
      attr.forEach((a) => {
        if (d[a] === undefined) throw "GroupedViolin: input data error.";
      });
      // if (0 == d.values.length) throw "Violin: Input data error";
    });
  }

  _addGroupDivider(dom) {
    const groups = this.scale.x.domain();
    const padding = Math.abs(
      this.scale.x(this.scale.x.domain()[1]) -
        this.scale.x(this.scale.x.domain()[0]) -
        this.scale.x.bandwidth()
    );

    const getX = (g, i) => {
      if (i !== groups.length - 1) {
        return this.scale.x(g) + +this.scale.x.bandwidth() + padding / 2;
      } else {
        return 0;
      }
    };

    dom
      .selectAll(".vline")
      .data(groups)
      .enter()
      .append("line")
      .classed("vline", true)
      .attr("x1", getX)
      .attr("x2", getX)
      .attr("y1", this.scale.y.range()[0])
      .attr("y2", this.scale.y.range()[1])
      .style("stroke-width", (g, i) => (i != groups.length - 1 ? 1 : 0))
      .style("stroke", "rgb(86,98,107)")
      .style("opacity", 0.5);
  }

  _validVertices(vertices) {
    return vertices.length > 0 && vertices.every(([x, y]) => Number.isFinite(x) && Number.isFinite(y));
  }
}
