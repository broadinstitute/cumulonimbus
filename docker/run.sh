#!/usr/bin/env bash
sudo docker build ~/git/cumulonimbus/docker/$1 -t $1
sudo docker run -it $1 bash
