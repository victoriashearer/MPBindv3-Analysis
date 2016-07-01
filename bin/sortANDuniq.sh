#!/bin/sh
Infile=$1
Outfile=$2
#convert fastq to txt
awk 'NR % 4 != 0' $Infile > tmp1
awk 'NR % 3 != 0' tmp1 > tmp2
tail -n +2 tmp2 > tmp && mv tmp tmp2
awk 'NR % 2 != 0' tmp2 > $Outfile
rm tmp1
rm tmp2
#remove redundant reads
sort $Outfile | uniq > tmp ; mv tmp $Outfile
#convert txt to fastq
#efficient way to add header line
awk '1;!(NR%1){print "@hES";}' $Outfile > tmp ; mv tmp $Outfile
string="@hES"; awk -v n=1 -v s="$string" 'NR == n {print s} {print}' $Outfile > tmp ; mv tmp $Outfile

#sed -i '1i @hES' $Outfile
#efficient way to add new line every 2 lines
awk '1;!(NR%2){print "+";}' $Outfile > tmp ; mv tmp $Outfile
# efficient way to add new line every 3 lines
awk '1;!(NR%3){print "blablablablablabalabkanalyzebetterlikethisbalabala";}' $Outfile > tmp ; mv tmp $Outfile
# efficient way to delete last line
head -n -1 $Outfile > tmp ; mv tmp $Outfile
#sed -i '$ d' $Outfile