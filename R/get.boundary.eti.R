#'@title Obtain the boundary of toxicity, efficacy and immune response of gBOIN-ETI design
#' @param target  a vector contains 3 target values for toxcity, efficay and immune reponse
#' @param p.tox.low the highest toxicity probability that is deemed subtherapeutic such that dose escalation should be made. The default value is 0.6.
#' @param p.tox.high the lowest toxicity probability that is deemed overly toxic such that deescalation is required. The default value is 1.4.
#' @param p.eff the percentage of deviation from the target, e.g., 0.6 indicates 40% deviation from target 
#' @param p.imm the percentage of deviation from the target, e.g., 0.6 indicates 40% deviation from target 
#' @details This function is to obtain the boundary of toxicity, efficacy and immune response of gBOIN-ETI design
#' @return \code{get.boundary.eti()} returns the operating characteristics of the gBOIN-ETI design as a data frame,
#'         including: (1) the boundary value of toxicity (\code{toxicity}),
#'         (2) the boundary value of efficacy (\code{efficacy}),
#'         (3) the boundary value of immune response (\code{immune})
#'         
#' @export 
get.boundary.eti <- function(target, # a vector contains 3 target values for toxcity, efficay and immune reponse
                             p.tox.low=0.6,
                             p.tox.high=1.4,
                             p.eff=0.6,
                             p.imm=0.6){	 	

  ########################### end of subroutines  ########################################
  
  # set cutoff target for hypothese
  tox.low=target[1]*p.tox.low
  tox.high=target[1]*p.tox.high
  eff=p.eff*target[2]
  immune=p.imm*target[3]
 
  # boundary for quasi-Bernoulli toxicity
  lambda.t1  = log((1-tox.low)/(1-target[1]))/log(target[1]*(1-tox.low)/(tox.low*(1-target[1])));
  lambda.t2  = log((1-target[1])/(1-tox.high))/log(tox.high*(1-target[1])/(target[1]*(1-tox.high)));
  
  # boundary for quasi-Bernoulli efficacy
  lambda.e  = log((1-eff)/(1-target[2]))/log(target[2]*(1-eff)/(eff*(1-target[2])));# < target[2]
  
  # boundary for continuous immune reponse
  lambda.imm = (target[3]+immune)/2
  
  out=list(toxicity=c(lambda.t1,lambda.t2),
           efficacy = lambda.e,
           immune=lambda.imm);
  
  return(out)
} 

get.boundary.eti(target=c(0.6,0.4,4),
                 0.6,1.4,0.6,0.6)
