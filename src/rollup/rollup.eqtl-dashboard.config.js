import nodeResolve from 'rollup-plugin-node-resolve';
import uglify from 'rollup-plugin-uglify';
import {minify} from 'uglify-es';

const name = 'EqtlDashboard';
export default {
    input:'src/' + name + '.js',
    output: {
        file: 'assets/js/eqtl-dashboard.bundle.min.js',
        format: 'iife',
        name: name
    },
    // sourcemap: 'inline',
    // name: name,
    plugins: [
        nodeResolve({jsnext: true, main: true}),
        uglify({}, minify)
    ]
}
