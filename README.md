
# ![RatGTEx](/images/RatGTExPortal2.png)

This repository contains code for the [RatGTEx web portal](https://ratgtex.org), which provides gene expression and eQTL data for different rat tissues.

---

### Using Rollup to bundle the applications
Most of our tools are written in ES6, with the exception of the GTEx Gene-eQTL Visualizer, and we recommend using a module bundler such as Rollup to recompile code if needed.

#### Rollup installation
To install Rollup and required libraries, you may run ```npm install``` in the repo's root directory on your computer. This will install the libraries under a subdirectory: node_modules.

#### Rollup configuration
The rollup configuration files for each tool is located in the directory [rollup](/rollup). To recompile a tool (e.g. GTEx Expression Map): run the following command in your local repo's root directory:

```rollup -c rollup/rollup.expression-map.config.js```

This will recompile and generate a new bundled tool code in the directory `build/js/`.

To minify the bundled code, first set the environment variable NODE_ENV to "prod", for example in a Bash terminal, the command would be:
```export NODE_ENV="prod"```

Then run rollup to recompile the code. All the demos are using the minified code of the tools.

#### Dependencies
RatGTEx Visualizations is distributed, in part, under and subject to the provisions of licenses for:

D3.js (https://d3js.org/), Copyright (c) 2017 Mike Bostock. All rights reserved.
Licensed under the BSD 3-clause license (https://opensource.org/licenses/BSD-3-Clause); and

jQuery (https://jquery.com/), Copyright (c) 2018 The jQuery Foundation. All rights reserved.
Licensed under the MIT license (https://jquery.org/license/).

### Acknowledgements
Much of this code is adapted from the [GTEx portal visualizations](https://github.com/broadinstitute/gtex-viz).
