#!/bin/bash

# This script executes a series of commands, replacing placeholders with two provided parameters,
# which are now accepted as command-line flags.

# Initialize variables for parameters
PARAM_VALUE1=""
PARAM_VALUE2=""

# Parse command-line flags
# -s: Server/Cluster/Application name (parameter1)
# -p: Pod/Deployment name (parameter2)
while getopts "s:p:" opt; do
  case $opt in
    s)
      PARAM_VALUE1="$OPTARG"
      ;;
    p)
      PARAM_VALUE2="$OPTARG"
      ;;
    \?) # Handle invalid options
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :) # Handle options requiring an argument but none provided
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Shift off the options and their arguments, so "$@" refers to the remaining non-option arguments.
shift $((OPTIND-1))

# Check if both parameters have been provided via flags
if [ -z "$PARAM_VALUE1" ] || [ -z "$PARAM_VALUE2" ]; then
  echo "Usage: $0 -s <server_name> -p <pod_name>"
  echo "Example: $0 -s my-application -p my-deployment"
  exit 1
fi

echo "Running commands with:"
echo "  Server Name (for nc, gcloud): $PARAM_VALUE1"
echo "  Pod Name (for kubectl exec): $PARAM_VALUE2"
echo "----------------------------------------------------"

# Command 1: Netcat verbose check
echo "Executing: nc -v \"$PARAM_VALUE1-gke.k8s.dev.bb.schrodinger.com\" 443"
nc -v "$PARAM_VALUE1-gke.k8s.dev.bb.schrodinger.com" 443 &

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: Netcat command failed. Exiting."
  exit 1
fi

echo "----------------------------------------------------"

# Command 2: Get GKE cluster credentials
echo "Executing: gcloud container clusters get-credentials \"$PARAM_VALUE1\" --region us-central1 --project livedesign-qa"
gcloud container clusters get-credentials "$PARAM_VALUE1" --region us-central1 --project livedesign-qa

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: gcloud command failed. Exiting."
  exit 1
fi

echo "----------------------------------------------------"

# Command 3: Display current kubectl context
echo "Executing: kubectl config current-context"
kubectl config current-context

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: kubectl config command failed. Exiting."
  exit 1
fi

echo "----------------------------------------------------"

# Command 4: Execute bash in a Kubernetes pod
echo "Executing: kubectl exec -n livedesign -it deploy/\"$PARAM_VALUE2\" -- bash"
kubectl exec -n livedesign -it deploy/"$PARAM_VALUE2" -- bash

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: kubectl exec command failed."
  # Note: kubectl exec might return non-zero if the bash session exits,
  # so you might want to handle this error check differently based on your needs.
fi

echo "----------------------------------------------------"
echo "Script execution completed."
