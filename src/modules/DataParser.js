/**
 * Copyright © 2015 - 2018 The Broad Institute, Inc. All rights reserved.
 * Licensed under the BSD 3-clause license (https://github.com/broadinstitute/gtex-viz/blob/master/LICENSE.md)
 */
"use strict";
export function getGtexUrls() {
  const host = "/api/v4/";
  return {
    singleTissueEqtl: host + "singleTissueEqtl?geneId=",
    dyneqtl: host + "dyneqtl?",
    variantId: host + "variant?variantId=",
    exon: host + "exon?geneId=",
    geneExp: host + "geneExpression?geneId=",
    medGeneExp: host + "medianGeneExpression?",
    topInTissueFiltered:
      host + "topExpressedGene?filterMtGene=true&tissueSiteDetailId=",
    topInTissue: host + "topExpressedGene?tissueSiteDetailId=",
    geneId: host + "gene?geneId=",
    tissue: host + "tissueInfo",
  };
}

/**
 * parse GTEx dyneqtl json
 * @param data {JSON} from GTEx dyneqtl web service
 * @returns data {JSON} modified data
 * @private
 */
export function parseDynEqtl(json) {
  // check required json attributes
  [
    "data",
    "genotypes",
    "pValue",
    "pValueThreshold",
    "tissueSiteDetailId",
  ].forEach((d) => {
    if (!json.hasOwnProperty(d)) {
      console.error(json);
      throw "Parse Error: Required json attribute is missing: " + d;
    }
  });

  json.expression_values = json.data.map((d) => parseFloat(d));
  json.genotypes = json.genotypes.map((d) => parseFloat(d));

  json.homoRefExp = json.expression_values.filter((d, i) => {
    return json.genotypes[i] == 0;
  });
  json.homoAltExp = json.expression_values.filter((d, i) => {
    return json.genotypes[i] == 2;
  });
  json.heteroExp = json.expression_values.filter((d, i) => {
    return json.genotypes[i] == 1;
  });

  // generate genotype text labels
  // let ref = json.variantId.split(/_/)[2];
  // let alt = json.variantId.split(/_/)[3];
  let ref = json.ref;
  let alt = json.alt;
  json.het = ref + alt;
  json.ref = ref + ref;
  json.alt = alt + alt;

  return json;
}

/**
 * Parse the single tissue eqtls from GTEx web service
 * @param data {Json}
 * @param tissueSiteTable {Json} optional for mapping tissueSiteDetailId to tissueSiteDetail, a dictionary of tissueSite objects (with the attr tissueSiteDetail) indexed by tissueSiteDetailId, and
 * @returns {List} of eqtls with attributes required for GEV rendering
 */
export function parseSingleTissueEqtls(data, tissueSiteTable = undefined) {
  const attr = "singleTissueEqtl";
  if (!data.hasOwnProperty(attr))
    throw "Parsing Error: required attribute is not found: " + attr;
  ["variantId", "tissueSiteDetailId", "nes", "pValue"].forEach((k) => {
    if (!data[attr][0].hasOwnProperty(k))
      throw "Parsing Error: required attribute is missing: " + attr;
  });

  return data[attr].map((d) => {
    d.x = d.variantId;
    d.displayX = d.variantId;
    d.y = d.tissueSiteDetailId;
    if (tissueSiteTable)
      d.displayY = tissueSiteTable[d.tissueSiteDetailId].tissueSiteDetail;
    d.value = d.nes;
    d.displayValue = d.nes.toPrecision(3);
    d.r = -Math.log10(d.pValue); // set r to be the -log10(p-value)
    d.rDisplayValue = parseFloat(d.pValue.toExponential()).toPrecision(3);
    return d;
  });
}

/**
 * Parse the genes from GTEx web service
 * @param data {Json}
 * @returns {List} of genes
 */
export function parseGenes(data, single = false, geneId = null) {
  const attr = "gene";
  if (!data.hasOwnProperty(attr))
    throw "Parsing Error: attribute gene doesn't exist.";
  if (data.gene.length == 0) {
    alert("Gene not found");
    throw "Fatal Error: gene(s) not found";
  }
  if (single) {
    if (geneId === null)
      throw "Please provide a gene ID for search results validation";
    if (data.gene.length > 1) {
      // when a single gene ID has multiple matches
      let filtered = data.gene.filter((g) => {
        return (
          g.geneSymbolUpper == geneId.toUpperCase() ||
          g.geneId == geneId.toUpperCase()
        );
      });
      if (filtered.length > 1) {
        alert("Fatal Error: input gene ID is not unique.");
        throw "Fatal Error: input gene ID is not unique.";
        return;
      } else if (filtered.length == 0) {
        alert("No gene is found with " + geneId);
        throw "Fatal Error: gene not found";
      } else {
        data.gene = filtered;
      }
    }
    return data.gene[0];
  } else return data[attr];
}

/**
 * Parse the tissues
 * @param data {Json}
 * @returns {List} of tissues
 */
export function parseTissues(json) {
  const attr = "tissueInfo";
  if (!json.hasOwnProperty(attr))
    throw "Parsing Error: required json attr is missing: " + attr;
  const tissues = json[attr];

  // sanity check
  ["tissueSiteDetailId", "tissueSiteDetail", "colorHex"].forEach((d) => {
    if (!tissues[0].hasOwnProperty(d))
      throw "Parsing Error: required json attr is missing: " + d;
  });

  return tissues;
}

/**
 * Parse the tissues and return a lookup table indexed by tissueSiteDetailId
 * @param json from web service tissueSiteDetail
 * @returns {*}
 */
export function parseTissueDict(json) {
  const attr = "tissueInfo";
  if (!json.hasOwnProperty(attr))
    throw "Parsing Error: required json attr is missing: " + attr;
  const tissues = json[attr];
  // sanity check
  ["tissueSiteDetailId", "tissueSiteDetail", "colorHex"].forEach((d) => {
    if (!tissues[0].hasOwnProperty(d))
      throw "Parsing Error: required json attr is missing: " + d;
  });
  return tissues.reduce((arr, d) => {
    arr[d.tissueSiteDetailId] = d;
    return arr;
  }, {});
}

/**
 * Parse the tissues sample counts, GTEx release specific
 * @param json from web service tissueInfo
 */
export function parseTissueSampleCounts(json) {
  const attr = "tissueInfo";
  if (!json.hasOwnProperty(attr))
    throw "Parsing Error: required json attr is missing: " + attr;
  const tissues = json[attr];

  // check json structure
  const tissue = tissues[0];
  if (!tissue.hasOwnProperty("tissueSiteDetailId"))
    throw "Parsing Error: required attr is missing: tissueSiteDetailId";
  if (!tissue.hasOwnProperty("rnaSeqAndGenotypeSampleCount"))
    throw "Parsing Error: required attr is missing: rnaSeqAndGenotypeSampleCount";
  return tissues;
}

/**
 * Parse the tissue groups
 * @param data {Json}
 * @param forEqtl {Boolean} restrict to eqtl tissues
 * @returns {Dictionary} of lists of tissues indexed by the tissue group name
 */
export function parseTissueSites(data, forEqtl = false) {
  // the list of invalid eqtl tissues due to sample size < 70
  // a hard-coded list because the sample size is not easy to retrieve
  const invalidTissues = [];

  const attr = "tissueInfo";
  if (!data.hasOwnProperty(attr))
    throw "Parsing Error: required json attribute is missing: " + attr;
  let tissues = data[attr];
  ["tissueSite", "tissueSiteDetailId", "tissueSiteDetail"].forEach((d) => {
    if (!tissues[0].hasOwnProperty(d))
      throw `parseTissueSites attr error. ${d} is not found`;
  });
  tissues =
    forEqtl == false
      ? tissues
      : tissues.filter((d) => {
          return !invalidTissues.includes(d.tissueSiteDetailId);
        }); // an array of tissueSiteDetailId objects

  // build the tissueGroups lookup dictionary indexed by the tissue group name (i.e. the tissue main site name)
  let tissueGroups = tissues.reduce((arr, d) => {
    let groupName = d.tissueSite + "Group"; // Added suffix to avoid ID collisions -DM
    let site = {
      id: d.tissueSiteDetailId,
      name: d.tissueSiteDetail,
    };
    if (!arr.hasOwnProperty(groupName)) arr[groupName] = []; // initiate an array
    arr[groupName].push(site);
    return arr;
  }, {});

  // modify the tissue groups that have only a single site
  // by replacing the group's name with the single site's name -- resulting a better Alphabetical order of the tissue groups

  Object.keys(tissueGroups).forEach((d) => {
    if (tissueGroups[d].length == 1) {
      // a single-site group
      let site = tissueGroups[d][0]; // the single site
      delete tissueGroups[d]; // remove the old group in the dictionary
      tissueGroups[site.name] = [site]; // create a new group with the site's name
    }
  });
  return tissueGroups;
}

/**
 * parse transcript isoforms from the GTEx web service: 'reference/transcript?release=v7&gencode_id='
 * @param data {Json} from web service exon
 * returns a list of all Exons
 */
export function parseExonsToList(json) {
  const attr = "exon";
  if (!json.hasOwnProperty(attr))
    throw "Parsing Error: required json attribute is missing: exon";
  return json[attr];
}

/**
 * parse median gene expression
 * @param data {Json} with attr medianGeneExpression
 * @returns {*}
 */
export function parseMedianExpression(data) {
  const attr = "medianGeneExpression";
  if (!data.hasOwnProperty(attr))
    throw "Parsing Error: required json attribute is missing: " + attr;
  const adjust = 1;
  // parse GTEx median gene expression
  // error-checking the required attributes:
  if (data[attr].length == 0) throw "parseMedianExpression finds no data.";
  ["median", "tissueSiteDetailId", "geneId"].forEach((d) => {
    if (!data[attr][0].hasOwnProperty(d)) {
      console.error(data[attr][0]);
      throw `Parsing Error: required json attribute is missingp: ${d}`;
    }
  });
  let results = data[attr];
  results.forEach(function (d) {
    d.value = Number(d.median);
    d.x = d.tissueSiteDetailId;
    d.y = d.geneId;
    d.displayValue = Number(d.median);
    d.id = d.geneId;
  });

  return results;
}

/**
 * parse the expression data of a gene for a grouped violin plot
 * @param data {JSON} from GTEx gene expression web service
 * @param colors {Dictionary} the violin color for genes
 */
export function parseGeneExpressionForViolin(
  data,
  useLog = true,
  colors = undefined
) {
  const attr = "geneExpression";
  if (!data.hasOwnProperty(attr))
    throw "Parsing Error: required json attribute is missing: " + attr;
  data[attr].forEach((d) => {
    ["data", "tissueSiteDetailId", "geneSymbol", "geneId"].forEach((k) => {
      if (!d.hasOwnProperty(k)) {
        console.error(d);
        throw "Parsing Error: required json attribute is missing: " + k;
      }
    });
    d.values = useLog
      ? d.data.map((dd) => {
          return Math.log10(+dd + 1);
        })
      : d.data;
    d.group = d.tissueSiteDetailId;
    d.label = d.geneSymbol;
    d.color = colors === undefined ? "#90c1c1" : colors[d.geneId];
  });
  return data[attr];
}

/* parse the expression data of a gene for boxplot
 * @param data {JSON} from GTEx gene expression web service
 * @param tissues {Object} mapping of tissue ids to labels (tissue name)
 * @param colors {Object} mapping of tissue ids to boxplot colors
 */
export function parseGeneExpressionForBoxplot(
  data,
  tissues = undefined,
  colors = undefined
) {
  const attr = "geneExpression";

  if (!data.hasOwnProperty(attr))
    throw `Parsing error: required JSON attribute ${attr} missing.`;

  data[attr].forEach((d) => {
    ["data", "geneId", "geneSymbol", "tissueSiteDetailId"].forEach((k) => {
      if (!d.hasOwnProperty(k)) {
        console.error(d);
        throw `Parsing error: required JSON attribute ${k} is missing from a record.`;
      }
    });
    d.label =
      tissues === undefined
        ? d.tissueSiteDetailId
        : tissues[d.tissueSiteDetailId];
    d.color = colors === undefined ? "#4682b4" : colors[d.tissueSiteDetailId];
  });

  return data[attr];
}
