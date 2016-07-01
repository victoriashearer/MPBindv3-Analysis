import commands,sys,os
from scipy import stats
from scipy.stats import norm
scriptPath=os.path.abspath(os.path.dirname(__file__));  ## absolute script path for all the modules

#e.g.  SELEX_file='H1_R5_SELEX.sense_seq.unique_sort.6mer'
#e.g.  Control_file='H1_R5_Control.sense_seq.unique_sort.6mer'
#e.g.  Outfile_name='H1_R5_SELEX_Vs_Control_Fisher_1_sided.txt'

SELEX_file=sys.argv[1]
Control_file=sys.argv[2]
Outfile_name=sys.argv[3]

###################################################
def function_SELEX_Vs_Control_Fisher_1_sided(SELEX_file,Control_file,Outfile_name):
# Initialize lists containing each k-mer with its correpsonding count data for the final selection and control rounds
    w_motif_SELEX={}  
    w_motif_Control={}
# Open k-mer count file for final selection round
    infile=open(SELEX_file,'r')
# Read header and ignore
    header_temp=infile.readline()
# Read second line
    line=infile.readline()
# While there are lines in the file:
    while(line):
# Assign value in first column to Motif variable (k-mer substring)
        Motif=line.split()[0]
# Assign k-mer count values to this variable
        Substring_counts=line.split()[1]
# Calculate the count of all other k-mers by subtracting Motif count from total count
        Substring_others=str(int(line.split()[2])-int(line.split()[1]))
# Assign fraction of Motif count
        Substring_fraction=line.split()[3]
# Assign all of the values to the Motif key in the list
        w_motif_SELEX[Motif]=[Substring_counts,Substring_others,Substring_fraction]
# Read next line and repeat
        line=infile.readline()
# Close file
    infile.close()
# Do the same thing for the control rounds
    infile=open(Control_file,'r')
    header_temp=infile.readline()
    line=infile.readline()
    while(line):
        Motif=line.split()[0]
        Substring_counts=line.split()[1]
        Substring_others=str(int(line.split()[2])-int(line.split()[1]))
        Substring_fraction=line.split()[3]
        w_motif_Control[Motif]=[Substring_counts,Substring_others,Substring_fraction]
        line=infile.readline()
    infile.close()

# Start writing to output file
    outfile=open(Outfile_name,'w')
# Write header row
    header='Motif'+'\t'+'P_value.Substring[1-sided]'+'\t'+'\t'+'Z-Score.Substring'+'\t'+'\n'
    outfile.write(header)

# For each k-mer in the count data list for final selection round    
    for k in w_motif_SELEX.keys():
# Write the k-mer string in the first column
        outfile.write(k+'\t')
# Assign its corresponding count values from the lists for the final selection and control rounds to variables
        Substring_a1=w_motif_SELEX[k][0] #SELEX kmer count
        Substring_a2=w_motif_SELEX[k][1] #SELEX other count
        Substring_b1=w_motif_Control[k][0] #Control kmer count
        Substring_b2=w_motif_Control[k][1] #Control other count
#       Substring_P_value=r.fisher_test(r.cbind([int(Substring_a1),int(Substring_a2)],[int(Substring_b1),int(Substring_b2)]),alternative='greater')['p.value']
#       Substring_Z_score=r.FunctionZScore(float(Substring_P_value))

# Run Fishers exact test using the four values for the k-mer (i.e. count, totalcount - count for both SELEX and control final round) 
        oddsratio, pval = stats.fisher_exact([[Substring_a1, Substring_a2], [Substring_b1, Substring_b2]], alternative='greater')
# Transform each p-value into a Z-score 
        Substring_P_value = pval
        if{pval < 10**(-15)}:
            z = 8.0    
        elif{pval > (1-10**(-15))}:
            z = -8.0
        else:
            z = norm.ppf(1-pval, loc=0, scale=1)            
        Substring_Z_score = z        
# Write the p-values and z-scores corresponding to the count-by-substring and count-by-read data to output file

        outfile.write(str(Substring_P_value)+'\t'+str(Substring_Z_score)+'\t'+'\n')
    outfile.close()

#############################################################################################
function_SELEX_Vs_Control_Fisher_1_sided(SELEX_file,Control_file,Outfile_name)
