
#' Enumeration of All Possible Four Integral Two-armed Bernoulli-outcomes. 
#'
#' This Function is to enumerate all possible two-armed Bernoulli-outcomes with the form of 2X2 tables When N patients to be treated. 
#' @param N        the number of patients to be treated.
#' @param X0       the initial information with setting c(0,0,0,0)
#' @return X=Xsm   a list of length N+1, each of which is a list of matrices. The values of matrices are 0 throughout.
#' @return M       the list of the output of 2-part of composition of intergers (n=0, 1,...,N) 
#' @return Z       a list with representing separately counts of A&success, A&failure, B&succsess, and B&failure; and the length of each list equals to the number of all possible combinations of Bernoulli-outcomes when N patients to be treated. 
#' @return sk      to get the same attributes of Xsm.  
#' @details   Decision process for the first patient starts with information 0, where no previous patients treated and so only one possible Bernoulli-outcoms with X0=c(0,0,0,0). 
#' @details   xsimple() function used here to enumerate all 2-part composition of N patients to be treated, on a simplex lattice (2,n=N). And this function is required to attach the package "combinat". 
#' @details   ZAS, ZAF,ZBS,and ZBF separately indicate the counts of  A&Success,A&Failure,B&Success, B&Failure.
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @keywords NA
#' @export
#' @examples
#install.packages("combinat")
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' attributes(out)


serialTables <- function(N,X0=rep(0,4)){
  M <- list()
  for(i in 0:N){
    M[[i+1]] <- as.matrix(xsimplex(2,i)) 
  }
  Xsm <- ZAS <- ZAF <- ZBS <- ZBF <- list()
  ZAS[[1]] <- matrix(X0[1],1,1)
  ZAF[[1]] <- matrix(X0[2],1,1)
  ZBS[[1]] <- matrix(X0[3],1,1)
  ZBF[[1]] <- matrix(X0[4],1,1)
  Xsm[[1]] <- list(matrix(0,1,1))
  for(i in 1:N){
    n <- i+1
    Xsm[[n]] <- list()
    ZAS[[n]] <- ZAF[[n]] <- ZBS[[n]] <- ZBF[[n]] <- list()
    for(j in 1:n){
      Xsm[[n]][[j]] <- matrix(0,j,n-j+1)
      ZAS[[n]][[j]] <- matrix(rep((j-1):0,n-j+1),j,n-j+1)
      ZAF[[n]][[j]] <- (j-1) - ZAS[[n]][[j]]
    }
    for(j in 1:n){
      ZBS[[n]][[j]] <- t(ZAS[[n]][[n+1-j]])
      ZBF[[n]][[j]] <- t(ZAF[[n]][[n+1-j]])
    }
  }
  Z <- list(unlist(ZAS),unlist(ZAF),unlist(ZBS),unlist(ZBF)) 
  sk <- attr(unlist(as.relistable(Xsm),recursive=TRUE),"skeleton")
  return(list(X=Xsm,M=M,Z=Z,sk=sk))
}


#' Calculate Exact Probability of a 2X2 Table with Bernoulli-outcomes
#'
#' This function is to calculate the exact probability of every possible 2X2 table with two-armed Bernoulli-outcomes enumerated in the serialTables () when N patients to be treated.
#' @param s1     a list of probabilities to select A-arm for the all possible 2X2 tables based on E.st_utinity () or T.st_utinity ().
#' @param p      a vector of true success rates of two arms, a and b. 
#' @return       a list of exact probabilities per 2X2 table. 
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#'
#' #Probability to select the A-arm based on utinity function of E.st: E.st_utinity()
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' E.table.probability<-table.prob(relist( E.st_prob.A,out$sk),p=c(0.8,0.6))
#' E.table.probability.v<-unlist(E.table.probability)
#' E.table.probability.v
#'
#' #Probability to select the A-arm based on utinity function of T.st: T.st_utinity() 
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(0.8,0.6))
#' T.table.probability.v<-unlist(T.table.probability)
#' T.table.probability.v
table.prob <- function(s1,p){
  p.A <- p[1]
  p.B <- p[2]
  ret <- s1
  ret[[1]][[1]] <- 1
  N <- length(s1)-1
  for(i in 1:N){
    n <- i+1
    ret[[n]] <- lapply(ret[[n]],"*",0)
    for(j in 1:i){
      tmp.select <- s1[[i]][[j]]
      dims <- dim(tmp.select)
      # (1) A & Success
      J <- j+1;xspan <- 1:dims[1];yspan <- 1:dims[2];
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + tmp.select * ret[[i]][[j]] * p.A
      # (2) A & Failure
      J <- j+1;xspan <- 2:(dims[1]+1);yspan <- 1:dims[2];
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + tmp.select * ret[[i]][[j]] * (1-p.A)
      # (3) B & Success
      J <- j;xspan <- 1:dims[1];yspan <- 1:dims[2];
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + (1-tmp.select) * ret[[i]][[j]]* p.B
      # (4) B & Failure
      J <- j;xspan <- 1:dims[1];yspan <- 2:(dims[2]+1);
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + (1-tmp.select) * ret[[i]][[j]] * (1-p.B)
    }
  }
  ret
}





#' Calculate the Fraction of A-arm Selection of Two-armed Bernoulli-outcomes.
#' 
#' Calculate the fraction of  A_armed for every 2X2 table numerated in the serialTables () when N patients to be treated. 
#' @param x     a vector of 4 elements including A$Success, A$Failure, B$Success, and B$Failure.
#' @details  A and B are examples of two arms, which resulting in favorable outcomes (success) and unfavorable outcomes (failure).
#' @return      the fraction of A-arm selected
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#' AS<-13
#' #  the successful counts of A treatment
#' AF<-5
#' #  the failures of A treatment
#' BS<-2
#' # the successful counts of B treatment
#' BF<-1
#' #the failures of B treatment
#' N<- AS+AF+BS+BF
#' A.frac(x=c(AS,AF,BS,BF))
#' 
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#' A.fraction<-apply(ABSF,1,A.frac)
#' # the fraction of A-arm selected per table.
#' A.fraction

A.frac <- function(x){
  if(sum(x)==0){
    0.5 
  }else{
    A.cont <- x[1]+x[2]
    Total.cont<- sum(x)
    A.cont/Total.cont
  }
}



#' Average Fraction of Favorable Outcomes of Two-armed Bernoulli-outcomes.
#' 
#' Calculate the average success fraction of two arms for every 2X2 table numerated in the serialTables () when N patients to be treated. 
#' @param x      a vector of 4 elements inlcuding A$Success, A$Failure, B$Success, and B$Failure.
#' @details   A and B are examples of two arms, which resulting in favorable outcomes (success) and unfavorable outcomes (failures).
#' @return       average of success fraction
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#' AS<-13
#' #  the successful counts of A treatment
#' AF<-5
#' #  the failures of A treatment
#' BS<-2
#' #  the successful counts of B treatment
#' BF<-1
#' #  the failures of B treatment
#' N<- AS+AF+BS+BF
#' success.frac(x=c(AS,AF,BS,BF))
#'
#'
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#' success.fraction<-apply(ABSF,1,success.frac)
#' #  the average fraction of successes per table.
#' success.fraction

success.frac <- function(x){
  if(sum(x)==0){
    0.5
  }else{
    fl <- x[2]+x[4]
    suc <- x[1]+x[3]
    suc/(suc+fl)
  }
  
}




#' Utinity Function Based on Expected Decision Strategy (E.st).
#'
#' This utinity function returns the probability to select one out of two-armed bandits with Bernoulli-outcomes (favorable and unfavorable) for the next (n+1)-th individual. 
#' E.st considerates the making-decision in a population with the homogeneously expected decision attitudes, which actually calculates the expected value of beta posterior distribution. 
#' @param x    the informed informaiton of n individuals (n=0,1,2,...) with consisted of 4 integers including counts of A$success, A$failure, B$success and B$failure. 
#' @details    A and B are examples of two treatment arms in the context, resulting into a binary outcomes success or failure after one patient to be treated. The four integeral outcomes are written as a vector.  
#' @details    In terms of every patient with respecting themseleves'decision attitudes, each decision process selecting either A or B treatment is a randomly probabilistical process which depends the output returned by their selected utinity function.   
#' @return prob_A   the probability to select A-arm, which is consisted of values 1, 0.5 or 0. 
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#' AS<-13
#'# the successful counts of A treatment
#' AF<-5
#'#the failures of A treatment
#' BS<-2
#'#the successful counts of B treatment
#' BF<-1
#'#the failures of B treatment
#' N<- AS+AF+BS+BF
#'#the total number of patients treated (n)
#' E.st_utinity(x=c(AS,AF,BS,BF))
#'# Decision to how much probability to select arm A for the (n+1)-th patient with an expected attitude. 

E.st_utinity<- function(x){ 
  x. <- x+1
  a <- x.[1]/(x.[1]+x.[2])
  b <- x.[3]/(x.[3]+x.[4])
  return(prob_A=(sign(a-b)+1)/2)
}




#' Utinity Function Based on Targeting Decision Strategy (T.st).
#' 
#' T.st considerates the making-decision in a population with target (optimistic/pessimistic) attitudes, and calculates the probability of favorable rate more than a targeting value. 
#' The optimism hope their targeting values higher than the higher expected value of Beta posterior out of two. In contrast, the pessimism hope their targeting values less likely than the higher expected one.
#' @param x        the informed informaiton of n individuals (n=0,1,2,...) with consisted of 4 integers including A$success, A$failure, B$success and B$failure counts. 
#' @param w        a degree of attitude of an individual,and is to parameterize his or her targeting value.  
#' @param w > 0    on behalf of an optimistic individual.
#' @param w < 0    on behalf of a pessimistic individual.
#' @details  A and B are examples of two treatment arms in the context, resulting into a binary outcomes success or failure after one patient to be treated. The four integeral outcomes are written as a vector.  
#' @details  In terms of every patient with respecting themseleves'decision attitudes, each decision process select either A or B treatment is a randomly probabilistical process which depends the output returned by their selected utinity function.   
#' @details  w range from -1 to 1. 
#' @return prob_A   the probability to select arm-A, which is consisted of values 1, 0.5 or 0. 
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#' 
#' AS<-13
#'# the successful counts of A treatment
#' AF<-5
#'#the failures of A treatment
#' BS<-2
#'#the successful counts of B treatment
#' BF<-1
#'#the failures of B treatment
#' N<- AS+AF+BS+BF
#'#the total number of patients treated (n)
#' w=0.5
#'# the degree of the (n+1)-th individual's attitude who is one of the optimism. 
#' T.st_utinity(x=c(AS,AF,BS,BF),w=0.5)
#'# Decide to how much probability to select arm A for the (n+1)-th patient with an optimistic attitude. 

T.st_utinity <- function(x,w){ 
  x. <- x+1
  a.exp <- x.[1]/(x.[1]+x.[2])
  b.exp <- x.[3]/(x.[3]+x.[4])
  tmp<-max(a.exp,b.exp)
  if(w>0){
    target<-tmp+(1-tmp)*w
  } else {
    target<-tmp*(1+w)
  }
  a <- pbeta(target,x.[1],x.[2],lower.tail=FALSE)
  b <- pbeta(target,x.[3],x.[4],lower.tail=FALSE)
  return(prob_A=(sign(a-b)+1)/2)
}




#' Weighed Average Favorable Outcomes Rate per Individual
#' 
#' Weighted mean success rate per individual, named as overall success rate, since it calculates the sum of the multiplication of the probability and average success rate for all possible 2X2 tables when N patients. 

#' @param xv      a vector of average success fraction of two arms in all possible tables.
#' @param pv      a vector of occurence probabilities of all tables
#' @param sk      a list of list structure of matrices indicating all possible Bernoulli-outcomes of two arms.
#' @return        a vector of overall success rate of a series of patients n=0,1,2,...,N.
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#'# all possible 2X2 tables
#'
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' #Probability to select the A-arm based on utinity function of E.st.
#' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(0.8,0.6))
#' E.table.probability.v<-unlist(E.table.probability)
#' success.fraction<-apply(ABSF,1,success.frac)
#' success.fraction[1]<-0.5
#' OSR(xv=success.fraction,pv=E.table.probability.v,sk=out$sk)
#'
#' # Probability to select the A-arm based on utinity function of T.st.
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(0.8,0.6))
#' T.table.probability.v<-unlist(T.table.probability)
#' success.fraction<-apply(ABSF,1,success.frac)
#' success.fraction[1]<-0.5
#' OSR(xv=success.fraction,pv=T.table.probability.v,sk=out$sk) 
OSR <- function(xv,pv,sk){
  wv <- xv * pv
  wl <- relist(wv,sk)
  sapply(wl,function(x)sum(unlist(x)))
}


#' Overall Fraction of A-arm Selection  per Individual
#' 
#' Average the fraction of A-arm selection per individual.

#' @param fv     a vector of the A-arm selection fraction out of two-armed Bernoulli in all possible tables.
#' @param pv     a vector of occurence probabilities of all tables
#' @param sk     a list of list structure of matrices with indicating all possible Bernoulli-outcomes when N patients treated.
#' @return       a vector of mean of A-arm selection fraction of a series of patients n=0,1,2,...,N.
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#' #all possible 2X2 tables
#'
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' #Probability to select the A-arm based on utinity function of E.st.
#' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(0.8,0.6))
#' E.table.probability.v<-unlist(E.table.probability)
#' A.fraction<-apply(ABSF,1,A.frac)
#' A.fraction[1]<-0.5
#' O_A_frac(fv=A.fraction,pv=E.table.probability.v,sk=out$sk)
#'
#' #Probability to select the A-arm based on utinity function of T.st.
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(0.8,0.6))
#' T.table.probability.v<-unlist(T.table.probability)
#' A.fraction<-apply(ABSF,1,A.frac)
#' A.fraction[1]<-0.5
#' O_A_frac(fv=A.fraction,pv=T.table.probability.v,sk=out$sk) 
O_A_frac<- function(fv,pv,sk){
  wv <- fv * pv
  wl <- relist(wv,sk)
  sapply(wl,function(x)sum(unlist(x)))
}



#' The Model of Homogeneous E.st
#' 
#' Return the overall success rate (OSR) based on the model of homogeneous expected decision attitudes (E.st).
#' @param  N         the number of patients to be treated
#' @param  a,b       the true success rates of two arms A and B treatment
#' @return  E.st_OSR         the Overall ratio of favorable outcomes(OSRs) for homogeneous E.st population.
#' @return  E.st_A.frac      the mean of the fraction of the A-arm treatment selected for homogeneous E.st population.
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples 
#'library(combinat) 
#' N<-10
#' a<-0.8
#' b<-0.6
#' Hom_E.st(N=N,a=a,b=b)

 Hom_E.st<-function(a,b,N){
  out<-serialTables(N)
  AS <- out$Z[[1]]
  AF <- out$Z[[2]]
  BS <- out$Z[[3]]
  BF <- out$Z[[4]]
  ABSF <- cbind(AS,AF,BS,BF)
  NS<-AS+AF+BS+BF
  Success.frac <- AS+BS/NS
  Success.frac[1]<- 0.5
  A.frac<-(AS+AF)/NS
  A.frac[1]<-0.5
  
  mat<-ABSF+1
  a.exp <- mat[,1]/(mat[,1]+mat[,2])
  b.exp <-  mat[,3]/(mat[,3]+mat[,4])
  E.st_prob_A<-(sign(a.exp-b.exp)+1)/2
  
  E.table.probability<-table.prob(relist(E.st_prob_A,out$sk),p=c(a,b))
  E.table.probability.v<-unlist(E.table.probability)
  E.st_OSR<-OSR(xv=Success.frac,pv=E.table.probability.v,sk=out$sk)
  E.st_A.frac<-O_A_frac(fv=A.frac,pv=E.table.probability.v,sk=out$sk)
  return(list(E.st_OSR=E.st_OSR,E.st_A.frac=E.st_A.frac))
 }





 
 #' The Model of Homogeneous T.st
 #' 
 #' Return the overall success rate (OSR) based on the model of homogeneous targeting decision attitudes (T.st).
 #' @param  N         the number of patients to be treated
 #' @param  a,b       the true success rates of two arms A and B treatment
 #' @param  w         the degree of decision attitudes. 
 #' @return  T.st_OSR         the Overall ratio of favorable outcomes(OSRs) for homogeneous T.st population.
 #' @return  T.st_A.frac      the mean of the fraction of the A-arm treatment selected for homogeneous T.st population.
 #' @seealso  \code{\link[combinat]{xsimplex}}
 #' @export
 #' @examples 
 #'library(combinat) 
 #' N<-10
 #' a<-0.8
 #' b<-0.6
 #' w<-0.5
 #' Hom_T.st(N=N,a=a,b=b,w=w)

Hom_T.st<-function(a,b,N,w){
  out<-serialTables(N)
  AS <- out$Z[[1]]
  AF <- out$Z[[2]]
  BS <- out$Z[[3]]
  BF <- out$Z[[4]]
  ABSF <- cbind(AS,AF,BS,BF)
  NS<-AS+AF+BS+BF
  Success.frac <- AS+BS/NS
  Success.frac[1]<- 0.5
  A.frac<-(AS+AF)/NS
  A.frac[1]<-0.5
  
  mat<-ABSF+1
  a.exp <- mat[,1]/(mat[,1]+mat[,2])
  b.exp <-  mat[,3]/(mat[,3]+mat[,4])
  tmp<-pmax(a.exp,b.exp)
  if(w>0){
    target<-tmp+(1-tmp)*w
  } else {
    target<-tmp*(1+w)
  }
  pa <- pbeta(target,mat[,1],mat[,2],lower.tail=FALSE)
  pb <- pbeta(target,mat[,3],mat[,4],lower.tail=FALSE)
  T.st_prob_A <- (sign(pa-pb)+1)/2
  T.table.probability<-table.prob(relist(T.st_prob_A,out$sk),p=c(a,b))
  T.table.probability.v<-unlist(T.table.probability)
  T.st_OSR<-OSR(xv=Success.frac,pv=T.table.probability.v,sk=out$sk)
  T.st_A.frac<-O_A_frac(fv=A.frac,pv=T.table.probability.v,sk=out$sk)
  return(list(T.st_OSR=T.st_OSR,T.st_A.frac=T.st_A.frac))
}




#' The Model of Heterogeneity in Decision Attitudes
#' 
#' Return the overall success rate (OSR) based on the model of heterogeneity decision attitudes (optimistic/pessimistic).
#' @param  iter      the times of randomly generated from W distribution. 
#' @param  N         the number of patients to be treated
#' @param  a,b       the true success rates of two arms A and B treatment
#' @param  u,v       two shape parameters of the distributions of decision attitudes w.  
#' @return            the mean of iter times of OSR
#' @seealso  \code{\link[combinat]{xsimplex}}
#' @export
#' @examples 
#'library(combinat) 
#'iter<-50
#' N<-10
#' a<-0.8
#' b<-0.6
#' u<-5
#' v<-5
#' Het_T.st(iter=iter,N=N,a=a,b=b,u=u,v=v)
Het_T.st<-function(iter,N,a,b,u,v){
  out<-serialTables(N)
  AS <- out$Z[[1]]
  AF <- out$Z[[2]]
  BS <- out$Z[[3]]
  BF <- out$Z[[4]]
  ABSF <- cbind(AS,AF,BS,BF)
  FL <- ABSF[,2]+ABSF[,4]
  SUC <- ABSF[,1]+ABSF[,3]
  Success.frac <- SUC/(SUC+FL)
  Success.frac[1]<- 0.5
  legs<-sapply(lapply(out$sk,function(x)sapply(x,length)),sum)
  legs.cum<-cumsum(legs)
  
  mean.success.per.n.het<-matrix(0,iter,(N+1))
  prob_A<-list()
  prob_A[[1]] <- 0.5
  for(i in 1:iter){
    for(j in 1:N){
      #Number of tables corresponds to each individual
      mat<-ABSF[(legs.cum[j]+1):legs.cum[j+1],]
      w<--1+2*rbeta(1,u,v)
      mat. <- mat+1
      a.exp <- mat.[,1]/(mat.[,1]+mat.[,2])
      b.exp <- mat.[,3]/(mat.[,3]+mat.[,4])
      tmp<-pmax(a.exp,b.exp)
      if(w>0){
        target<-tmp+(1-tmp)*w
      } else {
        target<-tmp*(1+w)
      }
      pa <- pbeta(target,mat.[,1],mat.[,2],lower.tail=FALSE)
      pb <- pbeta(target,mat.[,3],mat.[,4],lower.tail=FALSE)
      prob_A[[j+1]] <- (sign(pa-pb )+1)/2
    }             
    ####calculate OSR
    table.prob.v<-table.prob(relist(unlist(prob_A),out$sk),p=c(a,b))
    table.prob.het.v<-unlist(table.prob.v)
    mean.success.per.n.het[i,]<-OSR(xv=Success.frac,pv=table.prob.het.v,sk=out$sk)
    
  }
  
  het.mean.osr<-apply(mean.success.per.n.het,2,mean)
  return(het.mean.osr)
}




