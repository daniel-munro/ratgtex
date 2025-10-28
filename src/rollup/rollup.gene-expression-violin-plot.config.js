import nodeResolve from "rollup-plugin-node-resolve";
import terser from "@rollup/plugin-terser";

const name = "GeneExpressionViolinPlot";
export default {
  input: "src/" + name + ".js",
  output: {
    file: "assets/js/gene-expression-violin-plot.bundle.min.js",
    format: "iife",
    name: name,
  },
  plugins: [nodeResolve({ browser: true, main: true }), terser()],
};
