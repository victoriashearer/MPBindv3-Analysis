import commands,sys

Train_file=sys.argv[1]
Aptamer_file=sys.argv[2]
Sort_logical=sys.argv[3] # 'TRUE' or 'FALSE'
Output_file=sys.argv[4]
# A function that takes in the Z-scores for the list of k-mers in an aptamer sequence and calculates a combined Z-score
def function_combined_Z_score(Z_score_list):
# for each Z-score (i.e. for each k-mer)
    Z_score_list=[float(i) for i in Z_score_list]
# calculate the number of Z-scores (i.e. the number of k-mers)
    k=float(len(Z_score_list))
# calculate a combined score
    combined_Z_score=sum(Z_score_list)/(k**0.5)
    return combined_Z_score
# A function that takes in a list of all possible k-mers and their Z-scores, and an aptamer sequence and returns a combined Z-score for the aptamer
def function_aptamer_scan_meta_Z_score(w_motif_Z_score,aptamer_seq):
    #sequence length
    motif_len=len(w_motif_Z_score.keys()[0])
    substring_Z_score_list=[]
    # for each k-mer in the aptamer sequence
    for i in range(len(aptamer_seq)-motif_len+1):
        # append the k-mer Z-score to a list
        substring=aptamer_seq[i:i+motif_len]
        substring_Z_score_list.append(float(w_motif_Z_score[substring]))
    # Calculate a combined Z-score for an aptamer using the list
    meta_Z_score=function_combined_Z_score(substring_Z_score_list)
    substring_Z_score_list=[str(round(i,2)) for i in substring_Z_score_list]
    return ','.join(substring_Z_score_list)+'\t'+str(round(meta_Z_score,2))

# read motif Z_scores
w_motif_Z1={}
w_motif_Z3={}
w_motif_Z_Combined={}

infile=open(Train_file,'r')
header_temp=infile.readline()

line=infile.readline()
while(line):
    motif=line.split()[0]
    Z1=line.split()[1]
    Z3=line.split()[2]
    Combined_Z_Score=line.split()[3]
    w_motif_Z1[motif]=Z1
    w_motif_Z3[motif]=Z3
    w_motif_Z_Combined[motif]=Combined_Z_Score
    line=infile.readline()

infile.close()

if(Sort_logical.upper()=='FALSE'):
    infile=open(Aptamer_file,'r')
    outfile=open(Output_file,'w')
    header='Aptamer.Seq'+'\t'+'Z1.Scan'+'\t'+'Z1.MetaScore'+'\t'+'\t'+'Z3.Scan'+'\t'+'Z3.MetaScore'+'\t'+'Z_Combined.Scan'+'\t'+'Z_Combined.MetaScore'+'\n'
    outfile.write(header)
    line=infile.readline()
    while(line):
        Seq=line.split()[0].upper()
        outfile.write(Seq+'\t')
        outfile.write(function_aptamer_scan_meta_Z_score(w_motif_Z1,Seq)+'\t')
        outfile.write(function_aptamer_scan_meta_Z_score(w_motif_Z3,Seq)+'\t')
        outfile.write(function_aptamer_scan_meta_Z_score(w_motif_Z_Combined,Seq)+'\n')
        line=infile.readline()
    infile.close()
    outfile.close()
elif(Sort_logical.upper()=='TRUE'):
    out_list=[]
    infile=open(Aptamer_file,'r')
    line=infile.readline()
    while(line):
        Seq=line.split()[0].upper()
        this_Z1=function_aptamer_scan_meta_Z_score(w_motif_Z1,Seq)
        this_Z3=function_aptamer_scan_meta_Z_score(w_motif_Z3,Seq)
        this_combined=function_aptamer_scan_meta_Z_score(w_motif_Z_Combined,Seq)
        this_score=function_aptamer_scan_meta_Z_score(w_motif_Z_Combined,Seq).split()[-1]
        out_list.append([float(this_score),Seq+'\t'+this_Z1+'\t'+this_Z3+'\t'+this_combined+'\n'])
        line=infile.readline()
    infile.close()
    out_list.sort()
    out_list.reverse()
    outfile=open(Output_file,'w')
    header='Aptamer.Seq'+'\t'+'Z1.Scan'+'\t'+'Z1.MetaScore'+'\t'+'Z3.Scan'+'\t'+'Z3.MetaScore'+'\t'+'Z_Combined.Scan'+'\t'+'Z_Combined.MetaScore'+'\n'
    outfile.write(header)
    for i,j in out_list:
        outfile.write(j)
    outfile.close()
else:
    print 'Error ...'
    print 'Check -sort option ...'
    sys.exit(1)


