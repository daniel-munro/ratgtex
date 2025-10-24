import nodeResolve from "rollup-plugin-node-resolve";
import terser from "@rollup/plugin-terser";

const name = "BatchGeneExpression";
export default {
  input: "src/" + name + ".js",
  output: {
    file: "assets/js/expression-map.bundle.min.js",
    format: "iife",
    name: name,
  },
  plugins: [nodeResolve({ browser: true, main: true }), terser()],
};
