# input a1,a2,b1,b2

User_parameter_vector=commandArgs(trailingOnly=TRUE)
a1=as.integer(User_parameter_vector[1])
a2=as.integer(User_parameter_vector[2])
b1=as.integer(User_parameter_vector[3])
b2=as.integer(User_parameter_vector[4])

P_value=fisher.test(cbind(c(a1,a2),c(b1,b2)),alternative=c('greater'))$p.value

cat(P_value)

