import test, { beforeEach, afterEach } from "node:test";
import assert from "node:assert/strict";
import { JSDOM } from "jsdom";
import jquery from "jquery";
import { build, phenotypeLabel } from "../src/XqtlDashboard.js";

let dom, dashboard;
const modes = ["expression", "splicing"];
const tissues = ["Brain", "Liver"];
const urls = { tissue: "/tissues", geneId: "/gene?geneId=", phenotypes: "/phenotypes?", genotypePhenotypes: "/plots?" };
const json = (data, status = 200) => new Response(JSON.stringify(data), { status });
const item = (modality, tissue, phenotypeId) => ({ modality, tissueSiteDetailId: tissue, phenotypeId });
const inventory = () => ({
  datasets: modes.flatMap((modality) => tissues.map((tissueSiteDetailId) => ({ modality, tissueSiteDetailId, status: "ready" }))),
  phenotypes: [item("expression", "Brain", "GeneA"), item("expression", "Liver", "GeneA"),
    item("splicing", "Brain", "GeneA__chr1_10_20_clu_1_+"), item("splicing", "Liver", "GeneA__chr1_30_40_clu_1_+")],
});
const plots = (params) => ({
  genotypePhenotypes: [{ phenotypeId: params.get("modality") === "expression" ? "GeneA"
    : params.get("tissueId") === "Brain" ? "GeneA__chr1_10_20_clu_1_+" : "GeneA__chr1_30_40_clu_1_+",
  data: [-2, -1, 0, 1, 2], pValue: .001, associationStatus: "ok" }],
  genotypes: [0, 0, 1, 1, 2], ref: "A", alt: "G", pValueThreshold: .01,
});
const settle = async () => { for (let i = 0; i < 8; i++) await new Promise((resolve) => setTimeout(resolve, 0)); };

beforeEach(() => {
  dom = new JSDOM(`<details open><form id="xqtlQuery">
    <input name="gene" value="GeneA"><input name="variant" value="chr1:10">
    ${modes.map((m) => `<input type="checkbox" name="modality" value="${m}" data-label="${m}" checked>`).join("")}
    <div id="xqtlTissues"></div><button type="submit">Search</button></form></details>
    <div id="xqtlMessage"></div><div id="xqtlProgress"></div><div id="xqtlSummary"></div><div id="xqtlResults"></div>`);
  globalThis.document = dom.window.document;
  globalThis.window = dom.window;
  globalThis.XMLSerializer = dom.window.XMLSerializer;
  globalThis.$ = jquery(dom.window);
  dom.window.SVGElement.prototype.getBBox = () => ({ width: 15, height: 10 });
});
afterEach(() => { dashboard?.destroy(); dashboard = undefined; dom.window.close(); delete globalThis.IntersectionObserver; });

function mockRequests(t, handler) {
  t.mock.method(globalThis, "fetch", async (url, options) => {
    const parsed = new URL(url, "https://example.org");
    if (parsed.pathname === "/tissues") return json({ tissueInfo: tissues.map((id) => ({ tissueSiteDetailId: id, tissueSiteDetail: id })) });
    if (parsed.pathname === "/gene") return json({ gene: [{ geneId: parsed.searchParams.get("geneId") }] });
    return handler(parsed, options);
  });
}

test("dashboard aligns the phenotype union, shows real sample counts, and separates missing measurements", async (t) => {
  mockRequests(t, (url) => json(url.pathname === "/phenotypes" ? inventory() : plots(url.searchParams)));
  dashboard = await build({}, urls);
  await dashboard.search();
  await settle();
  assert.equal(document.querySelectorAll(".xqtl-modality").length, 2);
  assert.equal(document.querySelectorAll("tbody tr:not(.xqtl-modality)").length, 3);
  assert.equal(document.querySelectorAll("svg").length, 4);
  assert.equal([...document.querySelectorAll(".xqtl-empty")].filter((n) => n.textContent === "Not measured").length, 2);
  assert.match(document.querySelector(".xqtl-count").textContent, /n = 2/);
  assert.match(document.querySelector(".xqtl-significant").textContent, /\*/);
  const ids = [...document.querySelectorAll("[id]")].map((n) => n.id);
  assert.equal(new Set(ids).size, ids.length);
  assert.equal(document.querySelector("details").open, false);
});

test("inventory is displayed before plot requests and intersection loads a dataset only once", async (t) => {
  let callback;
  const observed = [];
  globalThis.IntersectionObserver = class {
    constructor(cb) { callback = cb; }
    observe(node) { observed.push(node); }
    disconnect() {}
  };
  let calls = 0;
  mockRequests(t, (url) => {
    if (url.pathname === "/phenotypes") return json(inventory());
    calls++;
    return json(plots(url.searchParams));
  });
  dashboard = await build({}, urls);
  await dashboard.search();
  assert.equal(calls, 0);
  assert.equal(observed.length, 4);
  callback([{ target: observed[0], isIntersecting: true }]);
  callback([{ target: observed[0], isIntersecting: true }]);
  await settle();
  assert.equal(calls, 1);
  assert.equal(document.querySelectorAll("svg").length, 1);
});

test("failed requests have retry controls and do not become missing measurements", async (t) => {
  let failed = false;
  mockRequests(t, (url) => {
    if (url.pathname === "/phenotypes") {
      const data = inventory();
      data.datasets.find((d) => d.modality === "splicing" && d.tissueSiteDetailId === "Liver").status = "unavailable";
      data.phenotypes = data.phenotypes.filter((p) => !(p.modality === "splicing" && p.tissueSiteDetailId === "Liver"));
      return json(data);
    }
    if (!failed && url.searchParams.get("tissueId") === "Brain" && url.searchParams.get("modality") === "expression") {
      failed = true;
      return json({}, 503);
    }
    return json(plots(url.searchParams));
  });
  dashboard = await build({}, urls);
  await dashboard.search();
  await settle();
  assert.match(document.getElementById("xqtlResults").textContent, /Dataset unavailable/);
  assert.match(document.getElementById("xqtlResults").textContent, /Could not load plot data/);
  const retry = [...document.querySelectorAll("button")].find((n) => n.textContent === "Retry");
  retry.click();
  await settle();
  assert.doesNotMatch(document.getElementById("xqtlResults").textContent, /Could not load plot data/);
});

test("replacing a search aborts old requests and prevents stale plots from appearing", async (t) => {
  let release, oldSignal;
  mockRequests(t, (url, options) => {
    if (url.pathname === "/phenotypes") return json(inventory());
    if (!release) {
      oldSignal = options.signal;
      return new Promise((resolve) => { release = () => resolve(json(plots(url.searchParams))); });
    }
    return json(plots(url.searchParams));
  });
  dashboard = await build({}, urls);
  await dashboard.search();
  document.querySelector('[name="gene"]').value = "GeneB";
  await dashboard.search();
  release();
  await settle();
  assert.equal(oldSignal.aborted, true);
  assert.match(document.getElementById("xqtlSummary").textContent, /^GeneB/);
  assert.equal(document.querySelectorAll("svg").length, 4);
});

test("single-pair validation runs before data retrieval", async (t) => {
  let calls = 0;
  mockRequests(t, () => { calls++; return json({}); });
  dashboard = await build({}, urls);
  document.querySelector('[name="gene"]').value = "GeneA,GeneB";
  await dashboard.search();
  assert.equal(calls, 0);
  assert.match(document.getElementById("xqtlMessage").textContent, /one gene/);
  assert.equal(phenotypeLabel("GeneA__XM_123.1", "isoforms", "GeneA"), "Transcript XM_123.1");
});

test("plot requests are limited to three concurrent tissue–modality datasets", async (t) => {
  const releases = [];
  let started = 0;
  mockRequests(t, (url) => {
    if (url.pathname === "/phenotypes") return json(inventory());
    started++;
    return new Promise((resolve) => { releases.push(() => resolve(json(plots(url.searchParams)))); });
  });
  dashboard = await build({}, urls);
  await dashboard.search();
  assert.equal(started, 3);
  releases.shift()();
  await settle();
  assert.equal(started, 4);
  releases.splice(0).forEach((release) => release());
  await settle();
});

test("SVG exports retain phenotype identity, genotype counts and median styling", async (t) => {
  mockRequests(t, (url) => json(url.pathname === "/phenotypes" ? inventory() : plots(url.searchParams)));
  let saved;
  globalThis.saveAs = (blob, name) => { saved = { blob, name }; };
  dashboard = await build({}, urls);
  await dashboard.search();
  await settle();
  document.querySelector(".xqtl-download").click();
  assert.match(saved.name, /GeneA.*Brain\.svg$/);
  const xml = await saved.blob.text();
  assert.match(xml, /GeneA · chr1:10/);
  assert.match(xml, /n = 2/);
  assert.match(xml, /stroke: #f88/);
});
