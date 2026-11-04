# ![RatGTEx](/assets/images/RatGTExPortal.png)

This repository contains code for the [RatGTEx web portal](https://ratgtex.org), which provides gene expression and eQTL data for different rat tissues.

There are also repositories for the RatGTEx [data pipeline](https://github.com/daniel-munro/ratgtex-pipeline) and to [process those results](https://github.com/daniel-munro/ratgtex-server-data) for this site's API and visualizations.

---

### Using Jekyll to build the site

This site uses [Jekyll](https://jekyllrb.com/) to build the pages for each genome version.

### Using Rollup to bundle the visualization scripts

Use Node.js 18 or later. Install the locked dependencies, then build all four bundles:

```shell
npm ci
npm run build
```

The generated scripts live in `assets/js/` and are ignored by Git. Rebuild them after editing `src/` and before building or deploying the Jekyll site.

To build just one visualization:

```shell
npx rollup -c src/rollup/rollup.expression-map.config.js
```

The individual configurations share `src/rollup/shared.js`. Their entry points are:

| Entry point | Used for |
| --- | --- |
| `GeneExpressionViolinPlot.js` | Gene page expression plots |
| `BatchGeneExpression.js` | Expression heatmaps and top expressed genes |
| `XqtlDashboard.js` | Multimodal xQTL dashboard |
| `GeneEqtlVisualizer.js` | Gene eQTL bubble map and genotype plots |

Reusable chart classes, parsers, tooltips, and SVG export helpers live in `src/modules/`. The page templates initialize the bundles and provide jQuery, Bootstrap, and FileSaver where needed. Keep chart selections scoped to their own SVG, and use D3's `(event, datum)` event handler signature.

Run the visualization regression checks with:

```shell
npm test
```

These tests exercise rendering and interactions with representative data in jsdom. They cover tooltip isolation, cluster hover, heatmap redraws, small sample groups, sample counts, empty eQTL results, and SVG downloads. jsdom does not perform browser layout, so use the local site to review responsive layout and appearance when changing styles or dimensions.

### API docs

See [api/README.md](api/README.md) for server configuration, expression indexing, API conventions, tests, and suggested next-release HPC changes. The API reads existing release files, uses per-gene phenotype indexes for genotype association queries, and queries exon/independent-association records through a SQLite catalog. Indexes are generated from those files without regenerating the scientific data.

Render the API doc HTML page from openapi.yaml (and add frontmatter specifying page URL):

```shell
npx @redocly/cli build-docs openapi.yaml -o pages/about_api.html
python3 - <<'PYDOC'
from pathlib import Path
path = Path("pages/about_api.html")
path.write_text("---\npermalink: /about/api/\n---\n{% raw %}\n" + path.read_text() + "\n{% endraw %}\n")
PYDOC
```

### Dependencies

RatGTEx Visualizations is distributed, in part, under and subject to the provisions of licenses for:

D3.js (https://d3js.org/), Copyright (c) 2017 Mike Bostock. All rights reserved.
Licensed under the BSD 3-clause license (https://opensource.org/licenses/BSD-3-Clause); and

jQuery (https://jquery.com/), Copyright (c) 2018 The jQuery Foundation. All rights reserved.
Licensed under the MIT license (https://jquery.org/license/).

### Acknowledgements

Visualization code was adapted from the [GTEx portal visualizations](https://github.com/broadinstitute/gtex-viz).
