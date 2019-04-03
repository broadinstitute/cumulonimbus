#!/usr/bin/env bash
project="broad-gdr-dig-storage"
name=$1
tag=$2
#if [[ -z ${project} ]]; then
#    echo "Need to specify name of Google project."
#    exit
#fi
if [[ -z ${name} ]]; then
    echo "Need to specify name of docker project ($(ls -dm */))."
    exit
fi
if [[ -z ${tag} ]]; then
    tag=$(date +%y%m%d)
fi
full="cumulonimbus-${name}:${tag}"
echo "Using Google project ${project}, Docker project ${name}, full tag ${full}"
echo "Submitting"
gcloud builds submit --tag gcr.io/${project}/${full} ~/git/cumulonimbus/docker/${name}
echo "Done"
