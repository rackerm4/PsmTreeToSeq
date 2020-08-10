FROM python:3.8

# Make directory for app
WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy sources
COPY . .

# Set entrypoint to pass arguments
ENTRYPOINT ["python", "/app/scripts/main.py"]
