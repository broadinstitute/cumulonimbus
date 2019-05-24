#!/usr/bin/env bash

version=$1

if [[ -z ${version} ]]; then
    echo "Need to specify chowser version."
    exit
fi

scp $HOME/git/chowser/target/chowser_${version}_all.deb silver:/web/personal/oliverr/software/chowser
