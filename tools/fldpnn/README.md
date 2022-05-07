# flDPnn Predictor
This repository contains the source files and binaries of the standalone version of [flDPnn](http://biomine.cs.vcu.edu/servers/flDPnn/).

**NOTE:** It is highly recommended that instead of using this repository, you use the docker version provided [here](https://gitlab.com/sina.ghadermarzi/fldpnn_docker). This way you can avoid complications arising from incompatible versions of requirements.


## Requirements
The following items are required for fldpnn to run. The versions on which we have tested fldpnn are in parantheses.
1. Linux (Ubuntu x64 20.04.2)  
2. tcsh  shell (6.21.00-1)
3. Java Runtime Enviornment (openjdk 1.0.8)
3. Python 3 (3.8.5)
4. Python packages:
	1. plotly (4.14.3)
	2. scikit-learn (0.23.1)
	3. keras (2.4.3)
	4. tensorflow (2.4.1)
	5. pandas (1.2.2)


## Running flDPnn
This program is only compatible with linux based enviornments. Before running the program, make sure that you unzip the dowloaded package inside a linux enviornment as well. 
flDPnn takes protein sequence(s) in a fasta file and produces outputs (text and visualizations) in the same folder as the input fasta file. An example input fasta with the desired outputs is in the folder `example`. To run flDPnn on the example input, run the following command in the main flDPnn folder:
> `python run_flDPnn.py example/test.fasta`


## Updates
- **2021 December** - Updated residue numbers in the output visualization. Updated output csv format.
- **2021 March** - Initial release
