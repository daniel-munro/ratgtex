/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
import { select } from "d3-selection";
import "d3-transition";

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

  show(info, event) {
    if (this.verbose) console.log(info);
    this.edit(info);
    this.move(event?.pageX, event?.pageY);
    select("#" + this.id)
      .interrupt()
      .style("display", "inline")
      .transition()
      .duration(this.duration)
      .style("opacity", 1.0);
  }

  hide() {
    select("#" + this.id)
      .interrupt()
      .transition()
      .duration(this.duration)
      .style("opacity", 0.0)
      .on("end", function () {
        select(this).style("display", "none");
      });
    this.edit("");
  }

  move(x, y) {
    // Callers pass coordinates from the D3 event; retain the last position
    // when showing a tooltip programmatically.
    x = x ?? this.x ?? 100;
    y = y ?? this.y ?? 100;
    this.x = x;
    this.y = y;

    if (this.verbose) {
      console.log(x);
      console.log(y);
    }
    x = x + this.offsetX;
    y = y + this.offsetY < 0 ? 10 : y + this.offsetY;
    select("#" + this.id)
      .style("left", `${x}px`)
      .style("top", `${y}px`);
  }

  edit(info) {
    select("#" + this.id).html(info);
  }
}
