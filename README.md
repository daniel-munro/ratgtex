# ![RatGTEx](/assets/images/RatGTExPortal.png)

This repository contains code for the [RatGTEx web portal](https://ratgtex.org), which provides gene expression and eQTL data for different rat tissues.

There are also repositories for the RatGTEx [data pipeline](https://github.com/daniel-munro/ratgtex-pipeline) and to [process those results](https://github.com/daniel-munro/ratgtex-server-data) for this site's API and visualizations.

---

### Using Jekyll to build the site

This site uses [Jekyll](https://jekyllrb.com/) to build the pages for each genome version.

### Using Rollup to bundle the visualization scripts

#### Rollup installation

To install Rollup and required libraries, you may run `npm install` in this directory. This will install the libraries under a subdirectory: `node_modules`.

#### Rollup configuration

The rollup configuration files for each visualization is located in the directory [src/rollup](/src/rollup). To recompile the code for a particular visualization, run the following command in this directory:

```shell
rollup -c src/rollup/rollup.expression-map.config.js
```

This will recompile and generate a new bundled visualization script in the directory `assets/js/`.

### API docs

Render the API doc HTML page from openapi.yaml (and add frontmatter specifying page URL):

```shell
npx @redocly/cli build-docs openapi.yaml -o pages/about_api.html
sed -i '' '1s|^|---\npermalink: /about/api/\n---\n|' pages/about_api.html
```

### Dependencies

RatGTEx Visualizations is distributed, in part, under and subject to the provisions of licenses for:

D3.js (https://d3js.org/), Copyright (c) 2017 Mike Bostock. All rights reserved.
Licensed under the BSD 3-clause license (https://opensource.org/licenses/BSD-3-Clause); and

jQuery (https://jquery.com/), Copyright (c) 2018 The jQuery Foundation. All rights reserved.
Licensed under the MIT license (https://jquery.org/license/).

### Acknowledgements

Visualization code was adapted from the [GTEx portal visualizations](https://github.com/broadinstitute/gtex-viz).
