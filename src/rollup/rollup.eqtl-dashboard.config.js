import nodeResolve from "rollup-plugin-node-resolve";
import terser from "@rollup/plugin-terser";

const name = "EqtlDashboard";
export default {
  input: "src/" + name + ".js",
  output: {
    file: "assets/js/eqtl-dashboard.bundle.min.js",
    format: "iife",
    name: name,
  },
  plugins: [nodeResolve({ browser: true, main: true }), terser()],
};
