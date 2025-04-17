#!/bin/bash

# Directory containing the folders
DIR="/eos/uscms/store/user/pprova/T5gg_production2/"

# List the directories and extract their names
for folder in $(ls -d "$DIR"*/ | grep -o '[^/]*/*$'); do
  echo "${folder%/}"
done
