#'@title Obtain the operating characteristics of A generalized Bayesian optimal interval design for dose optimization in immunotherapy by simulating trials.
#' @param target  a vector contains 3 target values for toxcity, efficay and immune reponse
#' @param prob.tox a 4*ndose matrix containing the true toxicity probabilities of each grade x the investigational dose levels.
#' @param prob.eff a 4*ndose matrix containing the true efficacy probabilities of each grade (PD,SD,PR,CR) the investigational dose levels.
#' @param mean.imm a vector with mean of imm response for each dose levels.
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param startdose the starting dose level for the trial 
#' @param p.tox.low the highest toxicity probability that is deemed subtherapeutic such that dose escalation should be made. The default value is 0.6.
#' @param p.tox.high the lowest toxicity probability that is deemed overly toxic such that deescalation is required. The default value is 1.4.
#' @param p.eff the percentage of deviation from the target, e.g., 0.6 indicates 40% deviation from target 
#' @param p.imm the percentage of deviation from the target, e.g., 0.6 indicates 40% deviation from target 
#' @param cv the coefficient of variation. e.g.,A CV of 25% was assumed to generate the immune response value
#' @param cutoff.eli.tox the cutoff to eliminate an overly toxic dose for safety. We recommend the default value of (\code{cutoff.eli.tox=0.95}) for general use.
#' @param cutoff.eli.eff the cutoff to eliminate an insufficient dose for efficacy. We recommend the default value of (\code{cutoff.eli.eff=0.99}) for general use.
#' @param ntrial the total number of trials to be simulated.
#' @param seed set seed
#' @param wt1 weight vector for calculating ETS
#' @param wt2 weight vector for calculating EES
#' @details This function is to obtain the operating characteristics of a generalized Bayesian optimal interval design for dose optimization in immunotherapy by simulating trials.
#' @return \code{get.oc.gBOIN_ETI()} returns the operating characteristics of the gBOIN-ETI design as a data frame,
#'         including: (1) selection percentage at each dose level (\code{selpercent}),
#'         (2) the number of patients treated at each dose level (\code{nptsdose}),
#'         (3) the number of toxicities observed at each dose level (\code{ntoxdose}),
#'         (4) the average number of toxicities (\code{ntox}),
#'         (4) the average number of efficacy (\code{neff}),
#'         (5) the average number of patients (\code{totaln}),
#'         (6) the percentage of early stopping without selecting the MTD (\code{pctearlystop}).   
#'         
#' @export 
get.oc.gBOIN_ETI <- function(target, # a vector contains 3 target values for toxcity, efficay and immune reponse
                             prob.tox, # a 4*ndose matrix
                             prob.eff, # a 4*ndose matrix 
                             mean.imm, # a vector with mean of imm response for each dose
                             ncohort, 
                             cohortsize, 
                             startdose=1, 
                             p.tox.low=0.6,
                             p.tox.high=1.4,
                             p.eff=0.6,
                             p.imm=0.6,
                             cv=0.25,
                             cutoff.eli.tox=0.95,
                             cutoff.eli.eff=0.99, # 0.99
                             ntrial=1000,
                             seed,
                             wt1=c(0,0,1,1),
                             wt2=c(0,0,1,1)){	 	
  ### simple error checking
  # dose skipping is not allowed
  
  ####################################################################################################
  ## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
  
  ### nned to update later
  select.mtd <- function(target.t,eff, y.t, y.e, n){
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function (x, wt = rep(1, length(x))) 
    {
      n <- length(x)
      if (n <= 1) 
        return(x)
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) 
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    ## determine whether the dose has been eliminated during the trial
    ndose=length(n);
    elimi=rep(0, ndose);
    for(i in 1:ndose)
    {
      if(n[i]>2) {if(1-pbeta(target.t, y.t[i]+1, n[i]-y.t[i]+1)>cutoff.eli.tox) {elimi[i:ndose]=1; break;}} #tox elimi
      if(n[i]>2) {if(pbeta(eff, y.e[i]+1, n[i]-y.e[i]+1)>cutoff.eli.eff) {elimi[i]=1;}} # eff elimi #lambda.e elimi[i:ndose]=1;
    }
    
    if(all(elimi[1:ndose]==1)) { selectdose=99; 
    } else{
      nadmis = which(elimi==0&n!=0); ## the highest admissble (or un-eliminated) dose level
      if (length(nadmis)==0) {selectdose=99}
      else if(length(nadmis)==1){
        selectdose=nadmis}
      else{
        ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior 
        phat = (y.t[nadmis]+0.005)/(n[nadmis]+0.01); 
        phat.var = (y.t[nadmis]+0.005)*(n[nadmis]-y.t[nadmis]+0.005)/((n[nadmis]+0.01)^2*(n[nadmis]+0.01+1))
        
        ## perform the isotonic transformation using PAVA
        phat = pava(phat, wt=1/phat.var) 
        phat = phat + (nadmis)*1E-10 ## break ties by adding an increasingly small number 
        selectdose = nadmis[sort(abs(phat-target.t), index.return=T)$ix[1]]  ## select dose closest to the target as the MTD
      }
    }
    return(selectdose);  
  }
  select.odimm <- function(target.t,eff,target.i, y.t, y.e,n,phat){
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function (x, wt = rep(1, length(x))) 
    {
      n <- length(x)
      if (n <= 1) 
        return(x)
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) 
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    ## determine whether the dose has been eliminated during the trial
    ndose=length(n);
    elimi=rep(0, ndose);
    for(i in 1:ndose)
    {
      if(n[i]>2) {if(1-pbeta(target.t, y.t[i]+1, n[i]-y.t[i]+1)>cutoff.eli.tox) {elimi[i:ndose]=1; break;}} #tox elimi
      if(n[i]>2) {if(pbeta(eff, y.e[i]+1, n[i]-y.e[i]+1)>cutoff.eli.eff) {elimi[i]=1;}} # eff elimi
    }
    #nadmis = which(phat>=target.i&elimi==0); #####????????????????????????????????????????????????????????????????????????????????
    nadmis = which(elimi==0&n!=0); #####large early termination
    
    ## perform the isotonic transformation using PAVA
    phat = pava(phat[nadmis], wt=rep(1, length(nadmis))) # 1 for all dose level weight
    phat = phat + (nadmis)*1E-10 ## break ties by adding an increasingly small number 
    selectdose = nadmis[sort(abs(phat-target.i), index.return=T)$ix[1]]  ## select dose closest to the target as the OD
    return(selectdose);  
  }
  ########################### end of subroutines  ########################################
  
  set.seed(seed);
  # set cutoff target for hypothese
  tox.low=target[1]*p.tox.low
  tox.high=target[1]*p.tox.high
  eff=p.eff*target[2]
  immune=p.imm*target[3]
  sigma.imm=mean.imm*cv
  
  ndose=ncol(prob.tox);	
  npts = ncohort*cohortsize;
  Y.tox=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  Y.eff=matrix(rep(0, ndose*ntrial), ncol=ndose) # store efficacy outcome
  Y.imm=matrix(rep(0, ndose*ntrial), ncol=ndose) # store imm response ###### alist
  
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  dselect = rep(0, ntrial); # store the selected dose level
  dselect.t = rep(0, ntrial); # store the selected dose level
  dselect.i = rep(0, ntrial); # store the selected dose level
  
  weight.T <- wt1
  weight.E <- wt2 # maybe modify the weight based on simulation setting.
  
  
  # boundary for quasi-Bernoulli toxicity
  lambda.t1  = log((1-tox.low)/(1-target[1]))/log(target[1]*(1-tox.low)/(tox.low*(1-target[1])));
  lambda.t2  = log((1-target[1])/(1-tox.high))/log(tox.high*(1-target[1])/(target[1]*(1-tox.high)));
  
  # boundary for quasi-Bernoulli efficacy
  lambda.e  = log((1-eff)/(1-target[2]))/log(target[2]*(1-eff)/(eff*(1-target[2])));# < target[2]
  
  # boundary for continuous immune reponse
  lambda.imm = (target[3]+immune)/2
  
  ################## simulate trials ###################
  for(trial in 1:ntrial)
  {
    tox.grade <- matrix(0,4,ndose) # grade of DLT (0,0.5,1,1.5)
    eff.grade <- matrix(0,4,ndose) # grade of eff (0,0.5,1,1.5)
    imm <- vector(mode = "list", length = ndose)
    
    y.t<-rep(0, ndose);    ## the number of DLT at each dose level
    y.e<-rep(0, ndose);    ## the number of eff at each dose level ???????????????????? y.imm
    y.i<-rep(0, ndose);
    #var.i<-rep(0, ndose);
    
    n<-rep(0, ndose);    ## the number of patients treated at each dose level
    mu.tox <- rep(0,ndose)
    mu.eff <- rep(0,ndose)
    mu.imm <- rep(0,ndose)
    
    earlystop=0;         ## indiate whether the trial terminates early
    d=startdose;         ## starting dose level
    elimi = rep(0, ndose);  ## indicate whether doses are eliminated
    
    for(i in 1:ncohort)  
    {  			
      ### generate grade toxicity outcome
      # leftp.tox <- 1-sum( prob.tox[,d]) # some prob left so that sum of prob is =1 in rmultinomial()
      # tox.grade[,d] <- tox.grade[,d]+rmultinom(1,cohortsize, prob=c(leftp.tox, prob.tox[,d]))[2:5]
      
      
      tox.grade[,d] <- tox.grade[,d]+rmultinom(1,cohortsize, prob=prob.tox[,d])
      eff.grade[,d] <- eff.grade[,d]+rmultinom(1, cohortsize, prob.eff[,d])
      imm[[d]] <-  c(imm[[d]], rnorm(n=cohortsize,mean=mean.imm[,d],sd=sigma.imm[,d]))
      
      n[d] = n[d] + cohortsize;
      #ETS
      y.t[d] = sum(tox.grade[,d]*weight.T)/max(weight.T)
      y.e[d] = sum(eff.grade[,d]*weight.E)/max(weight.E)
      y.i[d] = sum(imm[[d]])
      #var.i[d] =var(imm[[d]])
      
      mu.tox[d] <- y.t[d]/n[d]
      mu.eff[d] <- y.e[d]/n[d]
      mu.imm[d] <- y.i[d]/n[d]
      
      # determine which dose level should be eliminated based on tox
      
      if (1-pbeta(target[1], y.t[d]+1, n[d]-y.t[d]+1)>cutoff.eli.tox && n[d]>=3){ 
        elimi[d:ndose]=1;
        if (d == 1){
          earlystop = 1;
          break;
        }
      }
      
      
      if (pbeta(eff, y.e[d]+1, n[d]-y.e[d]+1)>cutoff.eli.eff && n[d]>=3){
        elimi[d]=1; # only this dose level is eliminated due to no efficacy observed, this early stop
      }
      
      
      
      adddose=NULL
      
      ## dose escalation/de-escalation
      if(mu.tox[d]<=lambda.t2){
        #1.1
        if (mu.eff[d]>=lambda.e & mu.imm[d]>=lambda.imm) dnext=d 
        #1.2
        if(mu.eff[d]>=lambda.e & mu.imm[d]<lambda.imm){
          if(d==ndose){
            dnext=d
          }else{
            adddose=c(d,d+1)
          }
        } 
        #1.3
        if(mu.eff[d]<lambda.e & mu.imm[d]>=lambda.imm){
          if(d==ndose){
            dnext=d
          }else{
            adddose=c(d,d+1)
          }
          
        }
        #1.4
        if(mu.eff[d]<lambda.e & mu.imm[d]<lambda.imm){
          if(d==ndose){
            dnext=d
          }else{
            dnext=d+1
          } 
        }
      }
      
      #3.1
      if(mu.tox[d]>lambda.t2){
        if(d==1){
          dnext=d 
        }else{
          dnext=d-1
        }
      } 
      
      
      if(length(adddose)!=0){
        if(n[max(adddose)]== 0){
          dnext=max(adddose)
        }else{
          pe <- y.e[adddose]/n[adddose]
          if(sum(pe==pe[which.max(pe)])==1){
            dnext=adddose[which.max(pe)]
          }else{
            mean.i=sapply(imm,sum)[adddose[pe==pe[which.max(pe)]]]/n[adddose[pe==pe[which.max(pe)]]]
            dnext=adddose[pe==pe[which.max(pe)]][which.max(mean.i)]
          }
        }
      }
      
      if (dnext>d & sum(elimi[(d+1):ndose]==0)!=0){ dnext=min(c((d+1):ndose)[elimi[(d+1):ndose]==0])
      } else if (dnext<d &sum(elimi[1:(d-1)]==0)!=0){dnext=max(c(1:(d-1))[elimi[1:(d-1)]==0])}
      d <- dnext
    }
    
    Y.tox[trial, ] = y.t
    Y.eff[trial, ] = y.e
    Y.imm[trial,] = mu.imm ################ 
    N[trial, ] = n
    if (earlystop == 1) {
      dselect.t[trial] = 99
    }else {
      dselect.t[trial]=select.mtd(target[1], eff,y.t,y.e, n);
    }
    
    dselect.i[trial]=select.odimm(target[1],eff,target[3],y.t,y.e,n, mu.imm); 
    
    
    #   
    if(!is.na(dselect.i[trial])& dselect.t[trial] !=99 && dselect.i[trial]<=dselect.t[trial]) {
      B=intersect(c(dselect.i[trial]:dselect.t[trial]), which(elimi==0)) 
    }else if(is.na(dselect.i[trial])){#|dselect.i[trial]>dselect.t[trial]) {   #####????????????????????????????????????????????????????????????????????????????????
      nadmis <- which(elimi==0&n!=0)
      B=nadmis[nadmis<=dselect.t[trial]]
    }else{ 
      B=NULL
    }
    if (length(B)==0){
      dselect[trial]=99
    }else{
      pe=y.e/n
      if (length(which(pe[B]==max(pe[B])))>1) {
        ind <- B[which(pe[B]==max(pe[B]))]
        mean.i=mu.imm[ind]/n[ind]
        dselect[trial]=ind[which(mean.i==max(mean.i))]        
      } else { 
        dselect[trial]= B[which.max(pe[B])] 
        
      }
    }
    
  }
  selpercent = rep(0, ndose)
  selpercent=rep(0, ndose);
  nptsdose = apply(N,2,mean);
  ntoxdose = apply(Y.tox,2,mean);
  neffdose= apply(Y.eff,2,mean);
  for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100; }
  
  out=list(selpercent=selpercent, patients=nptsdose, ntox=ntoxdose,neff=neffdose, totaltox=sum(Y.tox)/ntrial, totaln=sum(N)/ntrial,
           percentstop=sum(dselect== 99)/ntrial*100);
  
  return(out)
} 