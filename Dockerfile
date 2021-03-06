FROM debian:buster-slim
RUN mkdir -p /usr/share/man/man1 /usr/share/man/man2

RUN apt-get update && \
apt-get install -y --no-install-recommends make build-essential libssl-dev  wget curl llvm libidn11 openjdk-11-jre git nano tcsh sudo 

WORKDIR "/opt"

RUN git clone --depth=1 https://github.com/pyenv/pyenv.git .pyenv
ENV PYENV_ROOT="/opt/.pyenv" 
ENV PATH="${PYENV_ROOT}/shims:${PYENV_ROOT}/bin:${PATH}"

ENV PYTHON_VERSION=miniconda3-4.7.12
RUN pyenv install ${PYTHON_VERSION} && \
      pyenv global ${PYTHON_VERSION} && \      
      curl -sSL https://install.python-poetry.org | POETRY_HOME=/opt/poetry python3 - && \
      git clone https://github.com/wasicse/Dispredict3.0.git

WORKDIR "/opt/Dispredict3.0"
ENV PATH="/opt/poetry/bin:${PATH}"

ENV POETRY_VIRTUALENVS_IN_PROJECT="true"
RUN poetry install && \
      chmod -R 777 ../
ENTRYPOINT [ "/bin/bash" ]


