FROM phusion/baseimage:latest

# Make directory for app
WORKDIR /psm/

# Install dependencies
RUN apt-get update && apt-get install -y \
		git \
		seq-gen \
		python \
		python-pip \
		samtools \
	&& pip install dendropy \
	&& pip install -r requirements.txt \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Clone git
RUN git clone https://github.com/rackerm4/PsmTreeToSeq.git

# Create env
ENV PATH $PATH:/psm/PsmTreeToSeq/

ENTRYPOINT ["python", "/psm/main.py"]
