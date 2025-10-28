/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
import { select } from "d3-selection";

export default class Tooltip {
  constructor(
    id,
    verbose = false,
    offsetX = 30,
    offsetY = -40,
    duration = 100
  ) {
    this.id = id;
    this.verbose = verbose;
    this.offsetX = offsetX;
    this.offsetY = offsetY;
    this.duration = duration;
  }

  show(info) {
    if (this.verbose) console.log(info);
    this.edit(info);
    this.move();
    select("#" + this.id)
      .style("display", "inline")
      .transition()
      .duration(this.duration)
      .style("opacity", 1.0);
  }

  hide() {
    select("#" + this.id)
      .transition()
      .duration(this.duration)
      .style("opacity", 0.0);
    this.edit("");
  }

  move(x, y) {
    // If no coordinates provided, try to get them from the current event
    if (x === undefined || y === undefined) {
      // Try to get coordinates from the current mouse position
      if (typeof window !== "undefined" && window.event) {
        x = window.event.pageX;
        y = window.event.pageY;
      } else {
        // Fallback to default positioning
        x = 100;
        y = 100;
      }
    }

    if (this.verbose) {
      console.log(x);
      console.log(y);
    }
    x = x + this.offsetX; // TODO: get rid of the hard-coded adjustment
    y = y + this.offsetY < 0 ? 10 : y + this.offsetY;
    const t = select("#" + this.id)
      .style("left", `${x}px`)
      .style("top", `${y}px`);
  }

  edit(info) {
    select("#" + this.id).html(info);
  }
}
