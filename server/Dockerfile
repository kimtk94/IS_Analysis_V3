FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    python3 \
    python3-pip \
    python3-venv \
    r-base \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/is-analysis-venv
ENV PATH="/opt/is-analysis-venv/bin:${PATH}"

WORKDIR /app
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY . /app

CMD ["bash"]
