const path = require('path');

module.exports = {
    mode: 'production', // or 'development'
    entry: './src/webpack_main.js', // Your application's entry point
    output: {
        filename: 'node-bundle.js', // The name of your bundled file
        path: path.resolve(__dirname, 'lib'), // Output directory
        library: {
            name: 'pdd_node_bundle',
            type: 'var',
        },
    }
};