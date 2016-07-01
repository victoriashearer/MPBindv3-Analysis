import commands,sys

folder_path=sys.argv[1]  # e.g., /datadev/homes/pjiang/work/Aptamer_Pipeline/Peng_Test/Peng_PRAS_Out_U1/Unique_Reads/
RC=sys.argv[2]  # "NULL' or others
nmer=sys.argv[3] # 5,6

def function_combined_Z_score(Z_score_list):
    Z_score_list=[float(i) for i in Z_score_list]
    k=float(len(Z_score_list))
    combined_Z_score=sum(Z_score_list)/(k**0.5)
    return combined_Z_score

def function_combine_Z1_Z3(Z1_file,Z3_file,outfile):
    w_Z1={}  # Fisher.Substring
    w_Z3={}  # Spearman.Substring
    w_combine_Z={}
    # Z1
    infile=open(Z1_file,'r')
    header_temp=infile.readline()
    line=infile.readline()
    while(line):
        motif=line.split()[0]
        Z1=line.split()[-1]
        w_Z1[motif]=Z1
        w_combine_Z[motif]=''
        line=infile.readline()
    infile.close()
    # Z3
    infile=open(Z3_file,'r')
    header_temp=infile.readline()
    line=infile.readline()
    while(line):
        motif=line.split()[0]
        Z3=line.split()[-1]
        w_Z3[motif]=Z3
        line=infile.readline()
    infile.close()
    all_list=[]
    for i in w_combine_Z.keys():
        Z1=w_Z1[i]
        Z3=w_Z3[i]
        combined_Z_this=float(function_combined_Z_score([Z1,Z3]))
        out_line=i+'\t'+str(Z1)+'\t'+str(Z3)+'\t'+str(combined_Z_this)+'\n'
        all_list.append([combined_Z_this,out_line])
    all_list.sort()
    all_list.reverse()
    outfile=open(outfile,'w')
    header='Motif'+'\t'+'Z1[Fisher.Substring]'+'\t'+'Z3[Spearman.Substring]'+'\t'+'Combined_Z_Score'+'\n'
    outfile.write(header)
    for i,j in all_list:
        outfile.write(j)
    outfile.close()
if(RC!='NULL'):
    for i in nmer.split(','):
        Z1_file=folder_path+'Z1.'+i+'mer'
        Z3_file=folder_path+'Z3.'+i+'mer'
        outfile=folder_path+'Combined_Z_Score.train.'+i+'mer'
        function_combine_Z1_Z3(Z1_file,Z3_file,outfile)
else:
    for i in nmer.split(','):
        Z3_file=folder_path+'Z3.'+i+'mer'
        outfile=folder_path+'Combined_Z_Score.train.'+i+'mer'
        print Z3_file







