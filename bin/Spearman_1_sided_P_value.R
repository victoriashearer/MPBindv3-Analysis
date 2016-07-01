# input a1,a2

User_parameter_vector=commandArgs(trailingOnly=TRUE)
vector_1=as.numeric(strsplit(User_parameter_vector[1],split=',')[[1]])
vector_2=as.numeric(strsplit(User_parameter_vector[2],split=',')[[1]])

P_value=suppressWarnings(cor.test(vector_1,vector_2,method='spearman',alternative='greater')$p.value)

if(is.na(P_value))
  { P_value=1 } else {P_value=P_value}

cat(P_value)

