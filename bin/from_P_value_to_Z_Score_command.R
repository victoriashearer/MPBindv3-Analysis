# Usage: Rscript from_P_value_to_Z_Score_example.R <P-value>
# Out: Z-score

User_parameter_vector=commandArgs(trailingOnly = TRUE)
p_value=as.numeric(User_parameter_vector[1])

FunctionZScore <- function(P_value)
  {
    if(P_value<(10^(-15)))
        {Z_score=8.0}
    else if(P_value>(1-10^(-15)))
        {Z_score=(-8.0)}
    else
# Find the (1-pval)-th percentile of a normal distribution whose mean is zero and standard deviation is 1
        {Z_score=qnorm(p=(1-P_value),mean=0,sd=1,lower.tail=TRUE)} 
    Z_score
  }


Z_score=FunctionZScore(p_value)

cat (Z_score)


