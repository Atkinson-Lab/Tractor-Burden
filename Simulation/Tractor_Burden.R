args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}


Generate_dosage = function(sample, ancestry){
  return(sapply(sample, function(x){return(rbinom(1, size = x, prob = ancestry))}))
}


SIMULATOR <- function(param){ 
  num_sample = param[1]
  num_vars= param[2]
  admix_param1= param[3]
  admix_param2= param[4]
  AFR_maf= param[5]
  EUR_maf= param[6]
  num_risk_alleles= param[7]
  AFR_eff_size= param[8]
  EUR_eff_size= param[9]
  admix_eff_size= param[10]
  error_eff_size= param[11]
  
  admix = rbeta(num_sample, admix_param1, admix_param2)
  
  # AFR hapcount, generated from admix proportion with a binomial distribution (n = 2, prob = admix)
  AFR_hapcount = data.frame(t(sapply(admix, function(x){ return(rbinom(num_vars, size = 2, prob = x))})))
  EUR_hapcount = 2 - AFR_hapcount
  
  # generate AFR risk allele dosage
  AFR_dosage = data.frame(t(apply(AFR_hapcount,1, Generate_dosage, AFR_maf)))
  
  # generate EUR risk allele dosage
  EUR_dosage = data.frame(t(apply(EUR_hapcount,1, Generate_dosage, EUR_maf)))
  
  # generate TOT (total) copies of risk allele dosage
  TOT_dosage = AFR_dosage + EUR_dosage
  
  
  # determine which variants are risk alleles
  risk_alleles = sample(num_vars, size =num_risk_alleles, replace = F)
  admix_linear = log(admix/(1-admix))  # is this transformation necessary?
  error_term = rnorm(num_sample, sd = error_eff_size)
  
  
  # the outcome is a function of true risk alleles from EUR and AFR,
  # global ancestry (transformed into linear way), and an error term
  y = rowSums(AFR_dosage[,risk_alleles,drop = FALSE] * AFR_eff_size) + rowSums(EUR_dosage[,risk_alleles,drop = FALSE] * EUR_eff_size) + admix_linear * admix_eff_size + error_term
  
  
  # sum over AFR and EUR risk alleles, and perform a linear regression
  df = data.frame(cbind(rowSums(AFR_dosage), rowSums(EUR_dosage), rowSums(TOT_dosage), admix_linear, y))
  
  colnames(df) = c("AFR", "EUR", "TOT", "admix", "pheno")
  
  summary(lm(pheno ~ AFR + EUR + admix_linear, data = df))$coefficients[,4] 
  
  return(c(
    TOT = summary(lm(pheno~TOT + admix_linear, data = df))$coefficients[2,4],
    summary(lm(pheno ~ AFR + EUR + admix_linear, data = df))$coefficients[2:3,4]
  ))                     
}




params = crossing(num_sample = c(1000, 5000, 10000),
                  num_vars = 100, 
                  admix_param1 = 8, 
                  admix_param2 = 2,
                  AFR_maf = c(0.01, 0.05), 
                  EUR_maf = c(0.01, 0.05), 
                  num_risk_alleles = c(1, 5, 20),
                  AFR_eff_size = c(0, 1, 2), 
                  EUR_eff_size = c(0, 1, 2), 
                  admix_eff_size = 1,
                  error_eff_size = 2)



res = t(apply(params[1:10,], 1, SIMULATOR))

write.table(res, file, row.names = F, col.names = T, quote = F)




