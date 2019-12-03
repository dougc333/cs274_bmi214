#!/bin/bash

start=$SECONDS

python tanimoto.py data/drugs.csv data/targets.csv output_fast.csv

duration=$(( SECONDS - start ))

echo $duration
