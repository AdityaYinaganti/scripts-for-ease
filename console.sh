#!/bin/bash

CLUSTER_NAME="$1"

if [ -z "$CLUSTER_NAME" ]; then
  echo "Error: Please provide a cluster name as an argument."
  echo "Usage: ./$(basename "$0") <cluster_name>"
  exit 1
fi

# Function to open URL based on OS
open_url() {
  local url=$1
  if command -v xdg-open > /dev/null; then
    xdg-open "$url"
  elif command -v open > /dev/null; then
    open "$url"
  else
    echo "Please open your browser and go to: $url"
  fi
}

echo "Checking network connectivity to ${CLUSTER_NAME}-gke.k8s.dev.bb.schrodinger.com on port 443..."
nc -zv "${CLUSTER_NAME}-gke.k8s.dev.bb.schrodinger.com" 443

echo -e "\nGetting credentials for cluster '${CLUSTER_NAME}'..."
gcloud container clusters get-credentials "${CLUSTER_NAME}" --region us-central1 --project livedesign-qa

echo -e "\nChecking current kubectl context..."
kubectl config current-context

echo -e "\nOpening KOTS admin console..."
echo "The console usually runs at http://localhost:8800"

(sleep 3; open_url "http://localhost:8800") &

kubectl kots admin-console -n replicated
