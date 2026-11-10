# Introduction
See the program in action at [https://t-ext-explorer.com](https://t-ext-explorer.com).

This calculator computes the t-actions and maps in the exact sequence of Theorem 6.1 (see Section 4) of the paper `COMPUTING EXTENSIONS OF FLAG-JUMP t-MODULES`. (TODO: arXiv link).

See operator-focused README at [https://t-ext-explorer.com/t-Ext-Explorer-README.pdf](https://t-ext-explorer.com/t-Ext-Explorer-README.pdf) or in this repository at [https://github.com/dougkoellmer/t-ext-explorer/blob/main/t-Ext-Explorer-README.pdf](https://github.com/dougkoellmer/t-ext-explorer/blob/main/t-Ext-Explorer-README.pdf)

# Local Installation
- See `/scripts` folder for helpers.
- Run `cd /backend && source venv/bin/activate && pip install -r requirements.txt` to install Python dependencies.
- Use `scripts/local_api_server.sh` to start the local http server at `http://localhost:6789`.
- Use `scripts/local_server.sh` to start the local static file web server at `http://localhost:1234`.
- A local server can also be managed using Docker. See `/backend/Dockerfile`

# Local Installation Requirements
- Python 3 with pip.
- Any modern web browser.
- Node is used to build `/src/webpack_main.js` in case that need to be done.
- Only tested on Mac default Terminal application, but any POSIX-compliant shell should do.
- Jupyter notebook also provided in `/workbooks`.

# Architecture
- SymPy for all core math operations.
- Interactive HTML5/JavaScript front-end for ease of use.
- Front-end does some basic input validation, while back-end does the heavy lifting.
- Back-end can be run directly in various ways, without the front-end.
- Simple Python HTTP server with caching and multi-threading.
- Docker for deployment.
- Includes a Jupyter workbook that runs SymPy.

# Live Version Implementation Details
- API deployment handled through Docker.
- 60 second limit for calculations in live version hosted at [https://t-ext-explorer.com](https://t-ext-explorer.com).
- Live version is currently hosted on a rather slow server located in New York.
- To run larger examples more efficiently, it is recommended to install everything locally. See earlier instructions.
- Deployment of the HTML5 front-end is done through GitHub Pages using this very repository. Commits to the main branch get instantly reflected in the live version.

# Notes on Development History and Folder Organization
Development of this codebase took place over 9 months. For the first 6 months, the code was written manually, with all the mathematical operations done in such a way as to mimic how they would be done by hand. This early code played a role in shaping and refining the paper. However as the paper evolved, it eventually made the hand-written code (mostly) obsolete.

Once the paper approached a final version, we used ChatGPT to help generate the final reference implementation. See all files in `/backend/ai_code`. These ChatGPT-generated files were placed in a separate folder in order to make a clean distinction between the very different styles and conventions used in the original hand-written code.

The ChatGPT-generated code and the hand-written code were run side-by-side against thousands of automatically randomly generated examples to verify that both implementations produce the same output for Ext^1.

The hand-written code is still included in case it's helpful as a learning resource, and for any ongoing verification purposes. See all python files starting with `tee_`.

The hand-written code still calculates Ext^1 correctly, but does not produce any other necessary output.

# Automated Testing
- This aspect of the codebase could stand to be more formalized. 