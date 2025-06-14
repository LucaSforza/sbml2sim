FROM debian

WORKDIR /app

RUN apt-get update && apt-get install -y \
    make \
    python3 python3-pip python3.11-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /venv-sbml2sim
RUN echo "source /venv-sbml2sim/bin/activate" >> ~/.bashrc
RUN /venv-sbml2sim/bin/pip3 install nevergrad matplotlib requests

COPY . .

COPY ./lib64/libroadrunner.so /lib

RUN make -j4
