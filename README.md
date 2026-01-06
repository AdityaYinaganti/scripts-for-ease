# Scripts for ease
This repository contain scripts which make life easier for the QA to generate data and getting into different pods.

## Python script
Python file is creates a CSV file which creates different scientifically relevant HELM sequences which the user can use in different situations

## Requirements
* **Python**: 3.x

## Usage

### Running the script
Run the script with optional flags for length (`-l`) and number of sequences (`-n`).

python helm_maker.py -l 15 -n 100 -o my_peptides.csv


## Pod running and entering
This bash script can be used to get into the pod by giving 2 inputs

## Usage

 ### Bash script
 Run the script with necessary flags for server ('-s') and the pod name ('-p')
 
 bash server_and_pod.sh -p ld-accessories -s qa-demo-26-1

 ## Console logging
 This bash script can be used to login to the replicated console whenever the user wants to redeploy the server. 

 ### Running the script
 Run the sript with no flags but the name of the server you want to access the replicated server for.
 
 bash console.sh qa-demo-26-1

 ## Usage
 ### Python script
 Run the script with necessary flag to generate scientifically relevant SMILES string. 

 ### Requirements
 Install rdkit in the system to run the code successfully. 

 pip install rdkit

 ### Running the script
 Run the script with necessary flags for number of rows ('-r' or '--rows'), output file name ('-o' or '--output') and and and an optional flag for maximum length of the SMILES string ('-s' or '--size'). 

 python chemistry_maker.py --rows 10 --output my_file.csv --size 10 #--size is optional

 OR 

 python chemistry_maker.py -r 10 -o my_file.csv -s 10  #-s is optional

## Usage 
### Python script
Similar code to that of chemistry_maker.py but instead of using only certain elements, it will use the whole periodic table. 

### Requirements 
Install rdkit in the system to run the code successfully. 

pip install rdkit

### Running the script 
python chemistry_maker_alternate.py -r 10 -s 100 -o my_smiles.csv

OR 

python chemistry_maker_alternate.py --rows --s 100 --output my_smiles.csv
