# Dispredict3.0
Dispredict3.0: Prediction of Intrinsically Disordered Proteins with Protein
Language Model

### Table of Content

- [Setup](#getting-started)
- [Dataset](#Dataset)
- [Run with local OS](#Run-with-local-OS)
- [Prerequisites](#Prerequisites)
- [Download and install code](#download-the-code)
- [Run with Docker](#Run-with-Docker)
- [Run with Singularity](#Run-with-Singularity)
- [References](#References) 

# Getting Started
 

These instructions will get you a copy of the project up and running on your local machine or docker container for disorder prediction. 

 ## Dataset
The dataset can be found in the dataset directory. The train, test, and validation set is collected from [1].

## Run with local OS
### Download the code


- Retrieve the code

```
git clone https://github.com/wasicse/Dispredict3.0.git

```

### Prerequisites

We have tested Dispredict3.0 on Ubuntu 20.04. You would need to install the following software before replicating this framework in your local or server machine. 

1. pyenv latest version

    curl https://pyenv.run | bash
    exec $SHEL
    For more details, visit: https://github.com/pyenv/pyenv-installer

1. Python version 3.7.4

    pyenv install miniconda3-4.7.12
    pyenv local miniconda3-4.7.12 
    For more details, visit: https://github.com/pyenv/pyenv

2. Poetry version 1.1.13

    curl -sSL https://install.python-poetry.org | python3 - --version 1.1.13
    For more details, visit: https://python-poetry.org/docs/

3. tcsh  shell version 6.21.00-1

    git clone https://github.com/tcsh-org/tcsh
    cd tcsh
    ./configure
    make
    For more details, visit: https://github.com/tcsh-org/tcsh

3. Setup fldpnn tool

    Copy the tcsh executable file into the corresponding directory.

        cp tcsh ../Dispredict3.0/tools/fldpnn/programs/fMoRFpred/
        cp tcsh ../Dispredict3.0/tools/fldpnn/programs/DisoRDPbind/psipred  

### Run Dispredict3.0

To run the program, first install all required libraries by running the following command:

```
cd Dispredict3.0
poetry install
poetry shell
```

Then execute the following command to run Dispredict3.0 from the script directory.

```
cd script
poetry run python Dispredict3.0.py -f "../example/sample.fasta" -o "../output/"
```

- The following instructions show how to run dispredict3.0 with docker.

### Run with Docker
- To run the Dispredict3.0 tool with docker, you can either build the docker image using dockerfile or pull the docker image from the registry.
#### Build Docker image 

```
docker build -t wasicse/dispredict3.0 https://github.com/wasicse/Dispredict3.0.git#main    
```
 #### (Alternatively) Pull image from Docker registry.

- Pull the image from the registry.
 ```
docker pull wasicse/dispredict3.0
```
#### Run Dispredict3.0 using Docker image
- Create the dispredict3.0 container and mount the current (Dispredict3.0) directory (downlaoded from GitHub) into the docker container.

```
docker run -ti --name dispredict3.0  wasicse/dispredict3.0:latest
```

- Then, run following python commands inside the docker container to have the disordered prediction.

```
export PATH="/opt/poetry/bin:${PATH}"
source /opt/Dispredict3.0/.venv/bin/activate
python /opt/Dispredict3.0/script/Dispredict3.0.py -f "/opt/Dispredict3.0/example/sample.fasta" -o "/opt/Dispredict3.0/output/"
```

- Check **output** folder for results. The output should be available only inside the docker container. 

- You can also copy the output to the host computer using the following command:

```
docker cp dispredict3.0:/opt/Dispredict3.0/output/ .
```
### Run with Singularity 

- You can also run using Singularity using the following command.

```
singularity pull dispredict3.sif docker://wasicse/dispredict3.0
singularity run --writable-tmpfs dispredict3.sif
```
- Then, run following python commands inside the Singularity container to have the disordered prediction.

```
export PATH="/opt/poetry/bin:${PATH}"
source /opt/Dispredict3.0/.venv/bin/activate
python /opt/Dispredict3.0/script/Dispredict3.0.py -f "/opt/Dispredict3.0/example/sample.fasta" -o "/opt/Dispredict3.0/output/"
```

- The **output** folder should contain the results. The output directory contains the disorder probabilities with labels for each residue in **sample_disPred.txt** file. The fully disorder prediction for each protein sequence is stored in **sample_fullydisPred.txt** file.

## Authors

Md Wasi Ul Kabir, Md Tamjidul Hoque. For any issue please contact: Md Tamjidul Hoque, thoque@uno.edu 

## References

1. Hu, Gang, Akila Katuwawala, Kui Wang, Zhonghua Wu, Sina Ghadermarzi, Jianzhao Gao, and Lukasz Kurgan. “FlDPnn: Accurate Intrinsic Disorder Prediction with Putative Propensities of Disorder Functions.” Nature Communications 12, no. 1 (December 2021): 4438. https://doi.org/10.1038/s41467-021-24773-7.




