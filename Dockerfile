FROM python:3.12-slim

# System deps for matplotlib (basic) + certificates
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

WORKDIR /app

# Copy project metadata first for better layer caching
COPY pyproject.toml /app/pyproject.toml
COPY README.md /app/README.md

# Create venv + install deps
RUN uv sync --no-dev

# Copy the rest
COPY src /app/src
COPY app.py /app/app.py

EXPOSE 8501

# Streamlit needs to bind to 0.0.0.0 in containers
CMD ["uv", "run", "streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501"]
