#!/bin/sh
Infile=$1
Outfile=$2
awk '1;!(NR%1){print "@hES";}' $Infile > $Outfile
sed -i '1i @hES' $Outfile
awk '1;!(NR%2){print "+";}' $Outfile > tmp ; mv tmp $Outfile
awk '1;!(NR%3){print "blablablablablabalabkanalyzebetterlikethisbalabala";}' $Outfile > tmp ; mv tmp $Outfile
sed -i '$ d' $Outfile
