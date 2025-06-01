FROM debian

WORKDIR /app

COPY . .

COPY ./lib64/libroadrunner.so /lib
RUN export LD_LIBRARY_PATH=/lib64/

RUN apt-get update && apt-get install -y \
    make \
    python3 python3-pip python3.11-venv python3-tk \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv venv-sbml2sim
RUN venv-sbml2sim/bin/pip3 install -r requirements.txt

RUN make lib
