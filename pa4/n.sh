#!/bin/bash

start=$SECONDS

python networkgen.py data/drugs.csv data/targets.csv data/protein_nodes.csv

duration=$(( SECONDS - start ))

echo $duration
