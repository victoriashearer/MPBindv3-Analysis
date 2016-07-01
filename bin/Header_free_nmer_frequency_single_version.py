# Usage: python nmer_frequency_single_version.py <nmer> <infile> <outfile>
# E.g, python nmer_frequency_single_version.py 6 Round_0_unambiguous_SELEX_Seq_uppercase.txt 1_frequency_6mer_Round_0_SELEX.txt

import sys

if(len(sys.argv)!=4):
    print 'Usage: python nmer_frequency_single_version.py <nmer> <infile> <outfile>'
    sys.exit(1)
else:
    pass
#Parse args into module (k, infile, outfile)
user_defined_nmer=int(sys.argv[1])
infile_SELEX_Seq=sys.argv[2]
outfile_freq=sys.argv[3]

######## Function for Generate nmers #############################
## A function that takes in value of k and alphabet set and enumerates all possible k-mers into a list
def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

keys_tuple_nmer=list(product('ATCG', repeat=user_defined_nmer))
## A function that parses all possible k-mers into an array
def function_parse_nmer(keys_tuple_nmer):
    nmer_list_this=[]
    for i in keys_tuple_nmer:
        nmer_list_this.append(''.join(i))
    return nmer_list_this
# Array containing all possible k-mers
nmer_list=function_parse_nmer(keys_tuple_nmer)  # ***
# integer containing magnitude of k
motif_len=int(user_defined_nmer)
# A function that checks all characters in a read belong to the alphabet set of DNA (ATGC)
def Function_check_ATCG(read,w_ATCG):
    for i in range(len(read)):
        if(read[i] in w_ATCG):
            pass
        else:
            print 'Error: Non-{ATCG} charactor found in the sequence'
            print read
            exit(1)

w_ATCG={}
w_ATCG['A']=''
w_ATCG['T']=''
w_ATCG['C']=''
w_ATCG['G']=''

# print 'motif_len: '+str(motif_len)

#############################################################
#A hash map that contains all k-mer keys for count
w_nmer_by_substring={}
for i in nmer_list:
    w_nmer_by_substring[i]=int(0)
# A hash map that contains all k-mer keys for count-by-read
w_nmer_by_seq={}
for i in nmer_list:
    w_nmer_by_seq[i]=int(0)

##########################################
# A function that takes in value of k and a sequence and finds all k-mers present within the sequence
def function_split_into_substring(substring_len,seq):
    substring_list=[]
    for i in range(len(seq)-int(substring_len)+1):
        substring_list.append(seq[i:(i+int(substring_len))])
    return substring_list
#########################################
sum_all_seq_counts=int(0)

infile=open(infile_SELEX_Seq,'r')
# header_temp=infile.readline()

line=infile.readline()
#for all reads
while(line):
#split read into individual characters
    seq=line.split()[0].upper()
#check characters belong to alphabet set
    Function_check_ATCG(seq,w_ATCG)
#find all k-mers in read
    substring_list_this_seq=function_split_into_substring(motif_len,seq)
#for each k-mer found in read
    for k in substring_list_this_seq:
#increment count 
        w_nmer_by_substring[k]=w_nmer_by_substring[k]+1

    w_unique={}
#for each unique k-mer found in read
    for kk in substring_list_this_seq:
        w_unique[kk]=''
    for m in w_unique.keys():
#increment count-by-read
        w_nmer_by_seq[m]=w_nmer_by_seq[m]+1
#calculate total reads
    sum_all_seq_counts=sum_all_seq_counts+1
    line=infile.readline()

infile.close()

#### Outfile ###############
sum_all_substring_counts=int(0)
#calculate total k-mer count
for i in nmer_list:
    sum_all_substring_counts=sum_all_substring_counts+w_nmer_by_substring[i]
#write results to ouput file
outfile=open(outfile_freq,'w')
header='Motif'+'\t'+'Counts_by_substring'+'\t'+'Total_counts_all_substring'+'\t'+'Fraction_by_substring'+'\t'+'Counts_by_reads'+'\t'+'Total_reads'+'\t'+'Fraction_by_reads'+'\n'

outfile.write(header)

for k in nmer_list:
    outfile.write(k+'\t')
    Counts_by_substring=w_nmer_by_substring[k]
#calculate (k-mer count)/(total count) and write to file
    Fraction_by_substring=float(Counts_by_substring)/float(sum_all_substring_counts)
    outfile.write(str(Counts_by_substring)+'\t'+str(sum_all_substring_counts)+'\t'+str(Fraction_by_substring)+'\t')
#write count-by-read to file
    Counts_by_reads=w_nmer_by_seq[k]
#calculate (k-mer count-by-read)/(total reads) and write to file 
    Fraction_by_reads=float(Counts_by_reads)/float(sum_all_seq_counts)
    outfile.write(str(Counts_by_reads)+'\t'+str(sum_all_seq_counts)+'\t'+str(Fraction_by_reads)+'\n')

outfile.close()













