import nodeResolve from "rollup-plugin-node-resolve";
import terser from "@rollup/plugin-terser";

export function visualizationConfig(name, bundle) {
  return {
    input: `src/${name}.js`,
    output: {
      file: `assets/js/${bundle}.bundle.min.js`,
      format: "iife",
      name,
    },
    plugins: [nodeResolve({ browser: true, main: true }), terser()],
    onwarn(warning, warn) {
      // A missing import would leave a runtime dependency on a global variable.
      if (warning.code === "UNRESOLVED_IMPORT") throw new Error(warning.message);
      warn(warning);
    },
  };
}
