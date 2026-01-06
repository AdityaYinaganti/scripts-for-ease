# Scripts for ease
This repository contain scripts which make life easier for the QA to generate data and getting into different pods.

## Python script
Python file is creates a CSV file which creates different HELM sequences which the user can use in different situations

## Requirements
* **Python**: 3.x

## Usage

### Python Script
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

 ### Bash script
 Run the sript with no flags but the name of the server you want to access the replicated server for.
 
 bash console.sh qa-demo-26-1

