#!/usr/bin/env bash
cmd=$1
name=$2
tag=$3
if [[ -z ${cmd} ]]; then
    echo "Need to specify command."
    exit
fi
if [[ -z ${name} ]]; then
    echo "Need to specify name of docker project."
    exit
fi
if [[ -z ${tag} ]]; then
    tag=$(date +%y%m%d)
fi
full="cumulonimbus-${name}:${tag}"
echo "Using full tag ${full}"
echo "Building"
sudo docker build ~/git/cumulonimbus/docker/${name} -t ${full}
if [[ "${cmd}" == "build" ]]; then
    echo "Done building."
else
    if [[ "${cmd}" == "shell" ]]; then
        echo "Starting shell"
        sudo docker run -it ${full} bash
        echo "Done with shell"
    else
        echo "Unknown command ${cmd}."
    fi
fi
