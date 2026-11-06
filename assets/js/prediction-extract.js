(() => {
  const form = document.getElementById("prediction-extract");
  const status = document.getElementById("prediction-extract-status");
  const button = form.querySelector("button");

  function identifiers(value) {
    return [...new Set(value.trim().split(/[\s,]+/).filter(Boolean))];
  }

  function download(result) {
    const lines = [
      ["gene_id", "tissue", "modality", "phenotype_id", ...result.rfids].join("\t"),
      ...result.predictions.map((row) =>
        [row.geneId, row.tissue, row.modality, row.phenotypeId, ...row.values].join("\t"),
      ),
    ];
    const blob = new Blob([lines.join("\n") + "\n"], {
      type: "text/tab-separated-values;charset=utf-8",
    });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "ratgtex-predictions-round11_2.tsv";
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    const fields = new FormData(form);
    const genes = identifiers(fields.get("geneId"));
    const rfids = identifiers(fields.get("rfid"));
    status.classList.remove("text-danger");
    if (!genes.length || genes.length > 5 || !rfids.length || rfids.length > 100) {
      status.textContent = "Enter 1–5 gene symbols and 1–100 RFIDs.";
      status.classList.add("text-danger");
      return;
    }
    const params = new URLSearchParams({ geneId: genes.join(","), rfid: rfids.join(",") });
    for (const name of ["tissueId", "modality"]) {
      if (fields.get(name)) params.set(name, fields.get(name));
    }
    button.disabled = true;
    status.textContent = "Looking up predictions…";
    try {
      const response = await fetch(`/api/v4/predictions?${params}`);
      if (!response.headers.get("content-type")?.includes("application/json")) {
        throw new Error("Prediction lookup is unavailable.");
      }
      const result = await response.json();
      if (!response.ok) throw new Error(result.error || "Prediction lookup failed.");

      const hasResults = result.predictions.length && result.unmatchedRfids.length < result.rfids.length;
      if (hasResults) download(result);
      const messages = [hasResults ? `Downloaded ${result.predictions.length} phenotype rows.` : "No matching predictions."];
      if (result.unknownGenes.length) messages.push(`Unknown genes: ${result.unknownGenes.join(", ")}.`);
      if (result.genesWithoutPredictions.length) {
        messages.push(`No predictions in available datasets for: ${result.genesWithoutPredictions.join(", ")}.`);
      }
      if (result.unmatchedRfids.length) messages.push(`Unmatched RFIDs: ${result.unmatchedRfids.join(", ")}.`);
      if (result.unavailableDatasets.length) {
        messages.push(`Some requested datasets are unavailable (${result.unavailableDatasets.length} tissue and modality combinations).`);
      }
      status.textContent = messages.join(" ");
    } catch (error) {
      status.textContent = error.message;
      status.classList.add("text-danger");
    } finally {
      button.disabled = false;
    }
  });
})();
