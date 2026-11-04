import { select } from "d3-selection";
import { median } from "d3-array";
import GroupedViolin from "./modules/GroupedViolin.js";
import { getGtexUrls, getGenotypePhenotypesUrl } from "./modules/DataParser.js";
import { downloadSvg } from "./modules/utils.js";

const COLORS = ["#487f9f", "#d0a14b", "#56a090"];
const STATUS = {
  outside_cis_window: "Outside the cis window",
  variant_below_mapping_filter: "Below the mapping allele-frequency or allele-count filter",
  no_called_genotypes: "No called genotypes",
  no_residual_genotype_variation: "No genotype variation after covariate adjustment",
  no_residual_phenotype_variation: "No phenotype variation after covariate adjustment",
  nonfinite_phenotype: "Missing phenotype values",
};
const element = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};
const key = (modality, tissue) => `${modality}/${tissue}`;

export function phenotypeLabel(id, modality, gene) {
  if (modality === "expression") return "Expression level";
  if (modality === "stability") return "RNA stability";
  const suffix = id.startsWith(`${gene}__`) ? id.slice(gene.length + 2) : id;
  if (modality === "isoforms") return `Transcript ${suffix}`;
  const intron = suffix.match(/^(chr[^_]+)_(\d+)_(\d+)_clu_(\d+)_([+-])$/);
  if (intron) return `Intron ${intron[1]}:${intron[2]}–${intron[3]} · cluster ${intron[4]} (${intron[5]})`;
  return suffix.replace(/^grp_(\d+)_(?:upstream|downstream)_/, "Group $1 · ");
}

async function request(url, signal) {
  const response = await fetch(url, { credentials: "include", signal });
  if (!response.ok) {
    const error = new Error(`Request failed (${response.status})`);
    error.status = response.status;
    throw error;
  }
  return response.json();
}

/** One pair, with a metadata inventory and bounded, abortable plot requests. */
export async function build({ formId = "xqtlQuery", rootId = "xqtlResults", tissueId = "xqtlTissues",
  messageId = "xqtlMessage", progressId = "xqtlProgress", summaryId = "xqtlSummary" } = {}, urls = getGtexUrls()) {
  const form = document.getElementById(formId);
  const root = document.getElementById(rootId);
  const messages = document.getElementById(messageId);
  const progress = document.getElementById(progressId);
  const summary = document.getElementById(summaryId);
  const tissueMenu = document.getElementById(tissueId);
  const submit = form.querySelector('button[type="submit"]');
  let generation = 0, abort, observer;
  let destroyed = false;
  const modalities = [...form.querySelectorAll('input[name="modality"]')].map((input) => ({ id: input.value, name: input.dataset.label }));
  const onSubmit = (event) => { event.preventDefault(); search(); };
  form.addEventListener("submit", onSubmit);
  submit.disabled = true;
  try {
    const { tissueInfo } = await request(urls.tissue);
    tissueMenu.replaceChildren();
    for (const tissue of tissueInfo) {
      const label = element("label", "form-check form-check-inline");
      const checkbox = element("input", "form-check-input");
      checkbox.type = "checkbox";
      checkbox.name = "tissue";
      checkbox.value = tissue.tissueSiteDetailId;
      checkbox.dataset.label = tissue.tissueSiteDetail;
      checkbox.checked = true;
      label.append(checkbox, element("span", "form-check-label", tissue.tissueSiteDetail));
      tissueMenu.append(label);
    }
    submit.disabled = false;
  } catch {
    messages.textContent = "Tissues could not be loaded. Reload the page to try again.";
  }

  async function search() {
    if (destroyed) return;
    const geneInput = form.elements.namedItem("gene").value.trim();
    const variant = form.elements.namedItem("variant").value.trim()
      .replace(/^chr/i, "chr").replace(/^chr([xym]):/i, (_, chromosome) => `chr${chromosome.toUpperCase()}:`);
    const tissues = [...form.querySelectorAll('input[name="tissue"]:checked')].map((input) => ({ id: input.value, name: input.dataset.label }));
    const selectedModes = modalities.filter((mode) => [...form.querySelectorAll('input[name="modality"]:checked')].some((input) => input.value === mode.id));
    if (!geneInput || /[,\s]/.test(geneInput) || !/^chr(?:\d+|X|Y|M):[1-9]\d*$/i.test(variant)) {
      messages.textContent = "Enter one gene and one variant in chromosome:position format (for example, chr5:153176222).";
      return;
    }
    if (!tissues.length || !selectedModes.length) {
      messages.textContent = "Select at least one tissue and one modality.";
      return;
    }
    const run = ++generation;
    abort?.abort();
    observer?.disconnect();
    abort = new AbortController();
    const signal = abort.signal;
    messages.textContent = "";
    root.replaceChildren();
    summary.textContent = "";
    progress.textContent = "Finding measured phenotypes…";
    root.setAttribute("aria-busy", "true");
    const current = () => run === generation && !destroyed;
    try {
      const lookup = await request(urls.geneId + encodeURIComponent(geneInput), signal);
      if (!current()) return;
      const gene = lookup.gene?.[0]?.geneId;
      if (!gene) throw new Error(`Gene not found: ${geneInput}`);
      const params = new URLSearchParams({ geneId: gene, tissueId: tissues.map((t) => t.id).join(",") });
      if (selectedModes.length === 1) params.set("modality", selectedModes[0].id);
      const inventory = await request(urls.phenotypes + params, signal);
      if (!current()) return;
      summary.textContent = `${gene} · ${variant}`;
      const measured = inventory.phenotypes.filter((p) => selectedModes.some((m) => m.id === p.modality));
      const units = new Map();
      const cells = new Set();
      let sequence = 0, running = 0, bound = 4;
      const queue = [];
      for (const mode of selectedModes) for (const tissue of tissues) {
        const dataset = inventory.datasets.find((d) => d.modality === mode.id && d.tissueSiteDetailId === tissue.id);
        units.set(key(mode.id, tissue.id), { mode, tissue, metadata: new Map(), cells: [],
          state: dataset?.status === "ready" ? "idle" : "unavailable" });
      }
      for (const item of measured) units.get(key(item.modality, item.tissueSiteDetailId))?.metadata.set(item.phenotypeId, item);
      const viewport = element("div", "xqtl-viewport");
      viewport.tabIndex = 0;
      viewport.setAttribute("role", "region");
      viewport.setAttribute("aria-label", "Phenotypes by tissue. Scroll horizontally for more tissues and vertically for more phenotypes.");
      const table = element("table", "xqtl-grid");
      const caption = element("caption", "visually-hidden", `${gene} and ${variant}: normalized phenotypes grouped by modality, with tissues in columns.`);
      const head = element("thead");
      const header = element("tr");
      const first = element("th", "xqtl-label", "Modality / phenotype");
      first.scope = "col";
      header.append(first);
      for (const tissue of tissues) {
        const th = element("th", "xqtl-tissue", tissue.name);
        th.scope = "col";
        header.append(th);
      }
      head.append(header);
      table.append(caption, head);
      viewport.append(table);
      root.append(viewport);
      const updateProgress = () => {
        if (!current()) return;
        const eligible = [...units.values()].filter((unit) => unit.metadata.size);
        const finished = eligible.filter((unit) => unit.state === "ready").length;
        const failed = [...units.values()].filter((unit) => ["error", "unavailable"].includes(unit.state)).length;
        const count = cells.size;
        progress.textContent = `${count} phenotypes · ${finished} of ${eligible.length} tissue–modality datasets loaded. Scroll through the grid to load more.${failed ? ` ${failed} datasets unavailable or failed.` : ""}`;
      };
      const paint = (cell) => {
        if (!current() || !cell.visible || cell.painted) return;
        const unit = cell.unit;
        if (unit.state !== "ready") return;
        const phenotype = unit.data.genotypePhenotypes.find((p) => p.phenotypeId === cell.phenotype);
        cell.node.replaceChildren();
        if (!phenotype) {
          cell.node.append(element("p", "xqtl-empty", "Not measured"));
          cell.painted = true;
          return;
        }
        const entries = [0, 1, 2].map((dosage) => ({ group: "genotypes", label: String(dosage), color: COLORS[dosage],
          values: phenotype.data.filter((value, i) => Number.isFinite(value) && unit.data.genotypes[i] === dosage) }));
        entries.forEach((entry) => { entry.size = entry.values.length; });
        const max = Math.max(4, ...entries.flatMap((entry) => entry.values.map((value) => Math.ceil(Math.abs(value)))));
        if (max > bound) {
          bound = max;
          for (const row of cells) for (const other of row.cells) {
            if (other !== cell && other.painted) { other.painted = false; paint(other); }
          }
        }
        const validP = Number.isFinite(phenotype.pValue) && phenotype.associationStatus === "ok";
        const significant = validP && Number.isFinite(unit.data.pValueThreshold) && phenotype.pValue < unit.data.pValueThreshold;
        const top = element("div", "xqtl-cell-top");
        top.append(element("span", significant ? "xqtl-p xqtl-significant" : "xqtl-p",
          validP ? `p = ${phenotype.pValue.toExponential(2)}${significant ? " *" : ""}` : "Not tested"));
        cell.node.append(top);
        const svgId = `${rootId}-plot-${run}-${sequence++}`;
        const svg = select(cell.node).append("svg").attr("id", svgId).attr("width", 214).attr("height", 164)
          .attr("role", "img").attr("aria-label", `${cell.phenotype}, ${unit.tissue.name}. ${top.textContent}. Normalized phenotype by reference, heterozygous and alternate genotype.`);
        svg.append("title").text(`${gene} · ${variant} · ${cell.phenotype} · ${unit.tissue.name}. ${top.textContent}`);
        const dom = svg.append("g").attr("transform", "translate(36,10)");
        const violin = new GroupedViolin(entries);
        violin.render(dom, 168, 108, .1, undefined, [-bound, bound], "", false, 0, false, 0, false, false, false, false, false, false, 10);
        dom.selectAll(".violin-axis-label, .violin-x-axis, .violin-median").remove();
        dom.selectAll(".violin-g").each(function (entry) {
          select(this).selectAll(".violin-points circle").style("fill", entry.color);
        });
        const labels = unit.data.ref.length < 4 && unit.data.alt.length < 4
          ? [`${unit.data.ref}/${unit.data.ref}`, `${unit.data.ref}/${unit.data.alt}`, `${unit.data.alt}/${unit.data.alt}`]
          : ["Ref/Ref", "Ref/Alt", "Alt/Alt"];
        entries.forEach((entry, i) => {
          const x = violin.scale.subx(entry.label) + violin.scale.subx.bandwidth() / 2;
          dom.append("text").attr("x", x).attr("y", 127).attr("text-anchor", "middle").attr("class", "xqtl-genotype").text(labels[i]);
          dom.append("text").attr("x", x).attr("y", 142).attr("text-anchor", "middle").attr("class", "xqtl-count").text(`n = ${entry.size}`);
          if (entry.size) dom.append("line").attr("x1", x - 9).attr("x2", x + 9)
            .attr("y1", violin.scale.y(median(entry.values))).attr("y2", violin.scale.y(median(entry.values)))
            .attr("class", "violin-median").style("stroke", "#f88").attr("stroke-width", 2);
        });
        // GroupedViolin's IDs are local to a chart; give repeated small multiples unique IDs.
        dom.selectAll("[id]").attr("id", function () { return `${svgId}-${this.id}`; });
        dom.selectAll(".violin-g").append("title").text((entry) => `Genotype ${labels[Number(entry.label)]}: n = ${entry.size}, median = ${median(entry.values)?.toFixed(3)}`);
        const download = element("button", "btn btn-sm btn-link xqtl-download", "SVG");
        download.type = "button";
        download.setAttribute("aria-label", `Download ${cell.phenotype} in ${unit.tissue.name} as SVG`);
        download.addEventListener("click", () => downloadSvg($(svg.node()), `${gene}-${cell.phenotype}-${unit.tissue.id}.svg`, "unused-clone"));
        top.append(download);
        if (!validP) cell.node.append(element("p", "xqtl-note", STATUS[phenotype.associationStatus] || "Association statistic unavailable"));
        else if (!Number.isFinite(unit.data.pValueThreshold)) cell.node.append(element("p", "xqtl-note", "Significance threshold unavailable"));
        cell.painted = true;
      };
      const failed = (unit, error) => {
        unit.state = "error";
        for (const cell of unit.cells) {
          if (!unit.metadata.has(cell.phenotype)) continue;
          cell.node.replaceChildren(element("p", "xqtl-empty", error.status === 404 ? "Variant unavailable in this cohort" : "Could not load plot data"));
          const retry = element("button", "btn btn-sm btn-outline-secondary", "Retry");
          retry.type = "button";
          retry.addEventListener("click", () => { unit.state = "idle"; enqueue(unit); });
          cell.node.append(retry);
        }
      };
      const pump = () => {
        if (!current()) return;
        while (running < 3 && queue.length) {
          const unit = queue.shift();
          unit.state = "loading";
          running++;
          for (const cell of unit.cells) if (unit.metadata.has(cell.phenotype)) cell.node.textContent = "Loading…";
          request(getGenotypePhenotypesUrl(gene, variant, unit.tissue.id, unit.mode.id, urls), signal)
            .then((data) => {
              if (!current()) return;
              unit.data = data;
              unit.state = "ready";
              for (const cell of unit.cells) {
                if (!cell.visible && unit.metadata.has(cell.phenotype)) cell.node.textContent = "Scroll to view";
                paint(cell);
              }
            }).catch((error) => { if (current() && error.name !== "AbortError") failed(unit, error); })
            .finally(() => { running--; if (current()) { updateProgress(); pump(); } });
        }
      };
      const enqueue = (unit) => {
        if (unit.state !== "idle") return;
        unit.state = "queued";
        queue.push(unit);
        pump();
      };
      if (typeof IntersectionObserver !== "undefined") observer = new IntersectionObserver((events) => {
        for (const event of events) {
          const cell = event.target.xqtlCell;
          cell.visible = event.isIntersecting;
          if (cell.visible) { enqueue(cell.unit); paint(cell); }
        }
      }, { rootMargin: "200px" });
      let rowCount = 0;
      for (const mode of selectedModes) {
        const body = element("tbody");
        const modeRow = element("tr", "xqtl-modality");
        const modeHeading = element("th");
        modeHeading.colSpan = tissues.length + 1;
        modeHeading.scope = "rowgroup";
        modeRow.append(modeHeading);
        body.append(modeRow);
        const ids = [...new Set(measured.filter((p) => p.modality === mode.id).map((p) => p.phenotypeId))].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
        modeHeading.textContent = `${mode.name} · ${ids.length} ${ids.length === 1 ? "phenotype" : "phenotypes"}`;
        if (!ids.length) {
          const tr = element("tr");
          const label = element("th", "xqtl-label", mode.name);
          label.scope = "row";
          tr.append(label);
          for (const tissue of tissues) tr.append(element("td", "xqtl-empty", units.get(key(mode.id, tissue.id)).state === "unavailable" ? "Dataset unavailable" : "No measured phenotypes"));
          body.append(tr);
        }
        for (const id of ids) {
          rowCount++;
          const tr = element("tr");
          const label = element("th", "xqtl-label");
          label.scope = "row";
          label.append(element("span", "xqtl-phenotype-name", phenotypeLabel(id, mode.id, gene)), element("code", "xqtl-phenotype-id", id));
          tr.append(label);
          const row = { cells: [] };
          cells.add(row);
          for (const tissue of tissues) {
            const unit = units.get(key(mode.id, tissue.id));
            const td = element("td", "xqtl-cell");
            const cell = { node: td, unit, phenotype: id, visible: !observer, painted: false };
            td.xqtlCell = cell;
            row.cells.push(cell);
            unit.cells.push(cell);
            tr.append(td);
            if (unit.state === "unavailable") td.append(element("p", "xqtl-empty", "Dataset unavailable"));
            else if (!unit.metadata.has(id)) td.append(element("p", "xqtl-empty", "Not measured"));
            else { td.textContent = "Scroll to load"; observer?.observe(td); }
          }
          body.append(tr);
        }
        table.append(body);
      }
      root.setAttribute("aria-busy", "false");
      form.closest("details")?.removeAttribute("open");
      if (!rowCount) messages.textContent = "No measured phenotypes were found for these selections. Dataset availability is shown below.";
      if (!observer) for (const unit of units.values()) if (unit.metadata.size) enqueue(unit);
      updateProgress();
    } catch (error) {
      if (!current() || error.name === "AbortError") return;
      messages.textContent = error.status === 404 ? "The gene could not be found." : error.message.startsWith("Gene not found") ? error.message : "The phenotype list could not be loaded. Try the search again.";
      progress.textContent = "";
      root.setAttribute("aria-busy", "false");
    }
  }
  return { search, destroy() { destroyed = true; generation++; abort?.abort(); observer?.disconnect(); form.removeEventListener("submit", onSubmit); } };
}
