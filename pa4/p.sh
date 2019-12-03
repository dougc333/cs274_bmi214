#!/bin/bash

start=$SECONDS

python pvalue.py -n 500 -r 214 data/drugs.csv data/targets.csv P21918 P18089

duration=$(( SECONDS - start ))

echo $duration
