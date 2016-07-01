#!/bin/sh
total=$1
Infile=$2
awk '{$(NF+1)="'$total'";}1' OFS="\t" $Infile > tmp ; mv tmp $Infile
awk '{ if ($3 != 0) $(NF+1)=$2 / $3; else $(NF+1)=1 }1' OFS="\t" $Infile > tmp; mv tmp $Infile
sed -i '1iMotif\tCountsbyString\tTotalCounts\tRelativeFrequency' $Infile