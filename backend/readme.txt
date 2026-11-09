source venv/bin/activate
pip install -r requirements.txt




docker buildx build --platform linux/amd64 -t t-module-cruncher .
docker run -p 8080:8080 -e PORT=8080 t-module-cruncher

docker login
docker tag t-module-cruncher aytwit/t-module-cruncher:latest

docker push aytwit/t-module-cruncher:latest