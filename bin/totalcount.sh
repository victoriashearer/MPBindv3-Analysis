#!/bin/sh
Infile=$1
cat $Infile | awk '{sum+=$2 ; } END{print sum}'