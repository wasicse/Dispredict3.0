# Dispredict3.0
Dispredict3.0: Prediction of Intrinsically Disordered Proteins with Protein
Language Model

### Table of Content

- [Setup](#getting-started)
- [Dataset](#Dataset)
- [Prerequisites](#Prerequisites)
- [Download and install code](#download-and-install-code)
- [Run with local OS](#Run-with-local-OS)
- [Run with Docker](#Run-with-Docker)
- [References](#References) 

# Getting Started
 

These instructions will get you a copy of the project up and running on your local machine or docker container for disorder prediction. 

 ## Dataset
The dataset can be found in the dataset directory. The train, test, and validation set is collected from [1].

## Prerequisites

You would need to install the following software before replicating this framework in your local or server machine.

 ```
Python version 3.7.4

Poetry version 1.1.13
 ```
- Install poetry by running the following command:
 ```
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
 ```
To configure your current shell run the following command.
```
source $HOME/.poetry/env
 ``` 
## Download and install code

- Retrieve the code

```
git clone https://github.com/wasicse/Dispredict3.0.git
```

### Run with local OS

To run the program, first install all required libraries by running the following command:

```
cd Dispredict3.0
poetry install
```

Then execute the following command to run Dispredict3.0 from the script directory.

```
cd script
poetry run python Dispredict3.0.py -f "../example/sample.fasta" -o "../example/sample.fasta"
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
docker run -ti --name dispredict3.0 \
        -v /$(pwd):/home/$USER/Dispredict3.0 \
        wasicse/dispredict3.0:latest
```

- Then, run following python commands inside the docker container to have the disordered prediction.

```
source $HOME/.poetry/env
poetry shell
cd Dispredict3.0/script
poetry run python Dispredict3.0.py -f "../example/sample.fasta"
```

- Finally, check **output** folder for results. The output should be available in both the host and docker container. The output directory contains the disorder probabilities with labels for each residue in **sample_disPred.txt** file. The fully disorder prediction for each protein sequence is stored in **sample_fullydisPred.txt** file.



## Authors

Md Wasi Ul Kabir, Md Tamjidul Hoque. For any issue please contact: Md Tamjidul Hoque, thoque@uno.edu 

## References

1. Hu, Gang, Akila Katuwawala, Kui Wang, Zhonghua Wu, Sina Ghadermarzi, Jianzhao Gao, and Lukasz Kurgan. “FlDPnn: Accurate Intrinsic Disorder Prediction with Putative Propensities of Disorder Functions.” Nature Communications 12, no. 1 (December 2021): 4438. https://doi.org/10.1038/s41467-021-24773-7.




