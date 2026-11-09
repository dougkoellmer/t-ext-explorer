FROM python:3.10.16
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . /app
EXPOSE 8080
CMD ["python", "api_server.py"]

