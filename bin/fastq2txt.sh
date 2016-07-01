#!/bin/sh
Infile=$1
Outfile=$2
awk 'NR % 4 != 0' $Infile > tmp1
awk 'NR % 3 != 0' tmp1 > tmp2
sed -i 1d tmp2
awk 'NR % 2 != 0' tmp2 > $Outfile
rm tmp1
rm tmp2
