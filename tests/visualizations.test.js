import test, { beforeEach, afterEach } from "node:test";
import assert from "node:assert/strict";
import { JSDOM } from "jsdom";
import jquery from "jquery";
import { select } from "d3-selection";
import Tooltip from "../src/modules/Tooltip.js";
import GroupedViolin from "../src/modules/GroupedViolin.js";
import Heatmap from "../src/modules/Heatmap.js";
import Dendrogram from "../src/modules/Dendrogram.js";
import DendroHeatmap from "../src/modules/DendroHeatmap.js";
import { launch } from "../src/GeneExpressionViolinPlot.js";
import { render as renderEqtlViolin } from "../src/EqtlViolinPlot.js";
import { getGtexUrls, getGenotypePhenotypesUrl, parseSingleTissueEqtls } from "../src/modules/DataParser.js";
import { downloadSvg } from "../src/modules/utils.js";

let dom;
beforeEach(() => {
  dom = new JSDOM('<div id="plot"></div><div id="other"></div><div id="tip"></div>');
  globalThis.window = dom.window;
  globalThis.document = dom.window.document;
  globalThis.XMLSerializer = dom.window.XMLSerializer;
  globalThis.$ = jquery(dom.window);
  // jsdom has no layout engine; only label positioning uses this measurement.
  dom.window.SVGElement.prototype.getBBox = () => ({ width: 20, height: 10 });
});
afterEach(() => {
  select(document.documentElement).selectAll("*").interrupt();
  dom.window.close();
});
const svg = (id = "plot") => select(`#${id}`).append("svg").append("g");
const hover = (node) => node.dispatchEvent(new window.MouseEvent("mouseover", { clientX: 200, clientY: 150, bubbles: true }));
const values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

test("tooltip uses explicit event coordinates and survives rapid hide/show", async () => {
  const tip = new Tooltip("tip", false, 30, -40, 0);
  tip.show("First", { pageX: 200, pageY: 150 });
  assert.equal(document.querySelector("#tip").style.left, "230px");
  assert.equal(document.querySelector("#tip").style.top, "110px");
  tip.hide();
  tip.show("Second", { pageX: 250, pageY: 10 });
  await new Promise((resolve) => setTimeout(resolve, 40));
  assert.equal(document.querySelector("#tip").style.opacity, "1");
  assert.equal(document.querySelector("#tip").style.top, "10px");
  assert.equal(document.querySelector("#tip").textContent, "Second");
});

test("sample count labels support one or many entries with repeated counts", () => {
  const g = svg();
  for (const count of [1, 2, 3, 4]) {
    g.selectAll("*").remove();
    const plot = new GroupedViolin(Array.from({ length: count }, (_, i) => ({ group: "Tissue", label: `group${i}`, size: 11, values: [...values] })));
    plot.render(g, 500, 250, .05, undefined, [], "TPM", true, 30, true, 0, false, false, false, true);
    assert.deepEqual(g.selectAll(".violin-size-axis .tick text").nodes().map((n) => n.textContent), Array(count).fill("(11)"));
    assert.equal(g.selectAll("path.violin").size(), count);
  }
});

test("reset and scale updates preserve the small-sample point threshold", () => {
  const plot = new GroupedViolin([{ group: "Tissue", label: "AA", values: [1, 2, 3] }]);
  const g = svg();
  plot.render(g, 500, 250, .05, undefined, [], "TPM", true, 30, true, 0, false, false, false, false, false, false, 10);
  plot.reset();
  assert.equal(g.selectAll(".violin-points circle").size(), 3);
  plot.updateYScale("Expression", [0, 10]);
  assert.deepEqual(plot.scale.y.domain(), [0, 10]);
  assert.equal(g.selectAll(".violin-points circle").size(), 3);
  plot.updateXScale(["Tissue"]);
  assert.equal(g.selectAll(".violin-points circle").size(), 3);
});

test("heatmap redraw retains cells by coordinates, updates sizes and removes labels", () => {
  const data = [{ x: "Brain", y: "Gene1", value: 1 }, { x: "Brain", y: "Gene2", value: 1 }];
  const plot = new Heatmap(data, false);
  const g = svg();
  plot.draw(g, { w: 200, h: 100 });
  const cells = g.selectAll(".exp-map-cell").nodes();
  plot.data = data.map((d) => ({ ...d, value: 2 }));
  plot.redraw(g, ["Brain"], ["Gene2", "Gene1"], { w: 400, h: 200 });
  assert.deepEqual(g.selectAll(".exp-map-cell").nodes(), cells);
  assert.equal(+cells[0].getAttribute("width"), plot.xScale.bandwidth());
  assert.equal(+cells[0].getAttribute("y"), plot.yScale("Gene1"));
  assert.deepEqual(g.selectAll(".exp-map-ylabel").nodes().map((n) => n.textContent), ["Gene2", "Gene1"]);
  select("#other").classed("highlighted", true);
  hover(cells[0]);
  assert.match(document.querySelector("#heatmapTooltip").textContent, /Gene1/);
  plot.cellMouseout(g);
  assert.ok(document.querySelector("#other").classList.contains("highlighted"));
  plot.data = [plot.data[0]];
  plot.redraw(g, ["Brain"], ["Gene1"]);
  assert.equal(g.selectAll(".exp-map-cell").size(), 1);
  assert.equal(g.selectAll(".exp-map-ylabel").size(), 1);
});

test("dendrogram node hover works without a global d3 object", () => {
  for (const orientation of ["h", "v"]) {
    const g = svg();
    new Dendrogram("(Brain:1,Liver:2);", orientation).draw(g, 200, 100);
    const node = g.select(".dendrogram-node").node();
    hover(node);
    assert.equal(node.getAttribute("r"), "3");
  }
});

test("cluster hover receives the bound hierarchy node under D3 v7", () => {
  const g = svg();
  g.append("text").attr("class", "exp-map-xlabel").datum("Brain");
  const tree = new Dendrogram("(Brain:1,Liver:2);", "v");
  const chart = Object.create(DendroHeatmap.prototype);
  chart.visualComponents = { svg: g };
  const treeG = chart._renderTree("column", tree, { id: "tree", x: 0, y: 0, w: 200, h: 100 });
  hover(treeG.select(".dendrogram-node").node());
  assert.ok(g.select(".exp-map-xlabel").classed("highlighted"));
});

test("empty eQTL results return an empty array; malformed later records are rejected", () => {
  assert.deepEqual(parseSingleTissueEqtls({ significantCisEqtls: [] }), []);
  const record = { variantId: "chr1:1", tissueSiteDetailId: "Brain", nes: .3, pValue: .001 };
  assert.equal(parseSingleTissueEqtls({ significantCisEqtls: [record] })[0].r, 3);
  assert.throws(() => parseSingleTissueEqtls({ significantCisEqtls: [record, {}] }), /variantId/);
});

test("SVG export preserves applicable root, child and media CSS without modifying the chart", async () => {
  document.head.innerHTML = '<style>svg {font: 12px sans-serif} .violin {fill: red} .unused {fill: blue} @media print { .violin {stroke: black} }</style>';
  const chart = svg().append("path").attr("class", "violin");
  const original = chart.node().ownerSVGElement;
  const before = original.outerHTML;
  let saved;
  globalThis.saveAs = (blob, name) => { saved = { blob, name }; };
  downloadSvg($(original), "test.svg", "unused-clone");
  assert.equal(original.outerHTML, before);
  assert.equal(saved.name, "test.svg");
  const xml = await saved.blob.text();
  const parsed = new window.DOMParser().parseFromString(xml, "image/svg+xml");
  assert.equal(parsed.querySelector("parsererror"), null);
  assert.match(parsed.querySelector("style").textContent, /svg.*font/);
  assert.match(parsed.querySelector("style").textContent, /@media print/);
  assert.doesNotMatch(parsed.querySelector("style").textContent, /unused/);
});

test("gene expression tooltips remain attached to their own chart", async (t) => {
  t.mock.method(globalThis, "fetch", async (url) => new Response(JSON.stringify(url === "/tissue" ? {
    tissueInfo: [{ tissueSiteDetailId: "Brain", tissueSiteDetail: "Brain", colorHex: "123456" }]
  } : { geneExpression: [{ geneId: String(url), tissueSiteDetailId: "Brain", subsetGroup: "allData", unit: "TPM", data: [...values] }] })));
  const urls = { tissue: "/tissue", geneExp: "/gene/" };
  await launch("plot", "tip1", "Gene1", "Expression", urls, undefined, { w: 500, h: 250 });
  await launch("other", "tip2", "Gene2", "Expression", urls, undefined, { w: 500, h: 250 });
  const first = document.querySelector("#plot .violin-g");
  hover(first);
  assert.match(document.querySelector("#tip1").textContent, /Brain.*n = 11.*Median TPM: 5.000/);
  assert.equal(document.querySelector("#tip2").textContent, "");
  hover(document.querySelector("#other .violin-g"));
  assert.match(document.querySelector("#tip2").textContent, /Median TPM/);
  assert.ok(document.querySelector("#plot-spinner"));
  assert.ok(document.querySelector("#other-spinner"));
});


test("expression request URLs distinguish eQTL normalization from gene-page TPM", () => {
  const eqtl = new URL(getGenotypePhenotypesUrl("GeneA", "chr1:10", "Brain"), "https://example.org");
  assert.equal(eqtl.pathname, "/api/v4/genotypePhenotypes");
  assert.equal(eqtl.searchParams.get("modality"), "expression");
  assert.equal(eqtl.searchParams.get("variantId"), "chr1:10");
  const associations = new URL(getGtexUrls().significantCisEqtls + "GeneA", "https://example.org");
  assert.equal(associations.pathname, "/api/v4/significantCisEqtls");
  const geneExpression = new URL(getGtexUrls().geneExp + "GeneA", "https://example.org");
  assert.equal(geneExpression.searchParams.get("expressionScale"), "tpm");
});

test("eQTL violin requests and plots inverse-normal transformed sample values", async (t) => {
  const normalized = [-2, -1, 0, 1, 2];
  t.mock.method(globalThis, "fetch", async (url) => {
    const request = new URL(url, "https://example.org");
    assert.equal(request.pathname, "/api/v4/genotypePhenotypes");
    assert.equal(request.searchParams.get("modality"), "expression");
    return new Response(JSON.stringify({
      geneId: "GeneA", genotypePhenotypes: [{ phenotypeId: "GeneA", data: normalized, pValue: .01 }], genotypes: [0, 0, 1, 1, 2], ref: "A", alt: "G",
      pValue: .01, pValueThreshold: .05, tissueSiteDetailId: "Brain",
    }));
  });
  const par = {
    id: "plot", width: 500, height: 300,
    marginTop: 50, marginRight: 50, marginBottom: 50, marginLeft: 50,
    yLabel: "Norm. Expression",
  };
  await renderEqtlViolin(par, "GeneA", "chr1:10", "Brain");
  assert.deepEqual(par.data.map((group) => group.values), [[-2, -1], [0, 1], [2]]);
  assert.deepEqual(select("#plot").selectAll(".violin-points circle").data(), normalized);
});
