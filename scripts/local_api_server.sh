cd ../backend

# This is how to run the dev server if you need to include any Sage code.
# nodemon --exec "sage -python" api_server.py 6789

# This was how to run the server before inclusion of CompatibleBasis.py, which required Sage.
source venv/bin/activate
nodemon api_server.py 6789

# Variations on running the server directly, without hot reload. This was used early in development, but
# requires re-running the command every time a file changes, if I recall correctly.
# python3 -m http.server 1234
# python3 api_server.py 6789