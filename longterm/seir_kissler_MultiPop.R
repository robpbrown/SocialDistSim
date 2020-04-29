# Three parameters that determine the beta parameter for the SEIR
# R0
# distancing: is there social distancing: 0/1
# policyparameters: how does distancing affect R0. 
# Default: policyparameters$distancingeffect = 0.6
getbeta = function (t, distancing, R0, policyparameters) { 
    R = R0
    if (distancing){
        R = R * (1 - policyparameters$distancingeffect)
    }
    # Use R0 = beta/gamma for SIR. Use gamma  = 1/5 based on average of 5 days spent in I
    return (R/5)
}

getbeta_MultiPop = function (t, distancing, R0, policyparameters) { 
  R = R0
  if (distancing){
    R = R * (1 - policyparameters$distancingeffect)
  }
  # Use R0 = beta/gamma for SIR. Use gamma  = 1/5 based on average of 5 days spent in I
  return (R/5)
}

#S : 1
#E : 2
#I : 3
#HH: 4
#HC: 5
#C : 6
#R : 7
# KEY FUNCTION
# Make changes here to change policy or transition rates
maketransmatrix = function (t, currentstate, distancing, R0, policyparameters){
    beta = getbeta (t, distancing, R0, policyparameters)
    transmatrix = matrix (c (1,2, beta, 
                             2,3, 1/5,
                             3,4, 0.0308/5,
                             3,5, 0.0132/5,
                             3,7, 0.956/5,
                             4,7, 1/6,
                             5,6, 1/8,
                             6,7, 1/10), ncol = 3, byrow = T)
    return (transmatrix)
}


#S : 1
#E : 2
#I : 3
#HH: 4
#HC: 5
#C : 6
#R : 7
# KEY FUNCTION
# Make changes here to change policy or transition rates
maketransmatrix_MultiPop = function (t, currentstate, distancing, R0, policyparameters){
  beta = getbeta_MultiPop (t, distancing, R0, policyparameters)
  transmatrix = matrix (c (1,2, beta, 
                           2,3, 1/5,
                           3,4, 0.0308/5,
                           3,5, 0.0132/5,
                           3,7, 0.956/5,
                           4,7, 1/6,
                           5,6, 1/8,
                           6,7, 1/10), ncol = 3, byrow = T)
  return (transmatrix)
}

# KEY FUNCTION
# Use this to set SD
# One policy:  distancing is on if I(t) goes above i1 and off if I(t) goes below i2. Do nothing otherwise
distancingpolicy = function ( t, currentstate, distancing, policyparameters) { 
    if (currentstate[3]  >= policyparameters$hi){
        return (1)
    } 
    if (currentstate[3] <= policyparameters$lo){
        return (0)
    }
    return (distancing)
}




makeQmatrix = function (policyparameters, startt, stopt,nt, deltat,npops){
  #print(policyparameters)
   Qmat=matrix(0,npops,nt+1)
   pat=rep(c(rep(1,policyparameters$Qdur),rep(0,policyparameters$NQdur)),npops+ceiling(nt/(policyparameters$NQdur+policyparameters$Qdur)))
  
  # length(pat) 
  # print(pat)
   
   for (ii in c(1:nrow(Qmat))){
   #  print(pat[(((ii-1)*policyparameters$stagger)+1):(((ii-1)*policyparameters$stagger)+nt+1)])
   #  print(length(pat[(((ii-1)*policyparameters$stagger)+1):(((ii-1)*policyparameters$stagger)+nt+1)]))
   #  print(dim(Qmat))
     Qmat[ii,]=pat[(((ii-1)*policyparameters$stagger)+1):(((ii-1)*policyparameters$stagger)+nt+1)]
   }
  return (Qmat)
}

  
# Main function
# Note: In this SEIR model, we model the fraction of individuals in each compartment.
# 
# Inputs: 
#   startstate : matrix of column number of states rows number of populations. Each entry is >=0 and <=1. All entries sum to 1. Default c(PP,0,0,0,0,0,0) for each population where PP is the Pop Percentage of the sub pop
#   startt : starttime (in days) 
#   stopt: stoptime  (in days)
#   deltat: unit of time to simulate for one step. Default: 1/10 for 1/10 of a day. 
#   R0: reproduction number
#   policyparameters: list of all the things that are relevant for the policy. e.g. how much does distancing reduce R0
#   debug: for verbose mode
# 
# Outputs: A list. 
#   states: matrix of time X # number of states. Use this to plot e.g. infected vs time

#startstate=matrix(c(.33,.33,.34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=3,ncol=7)
#startstate=matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=3,ncol=7)
#startstate=matrix(c(0.33,0.33,0.33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=3,ncol=7)

#startt=1
#stopt=20
#R0=2
#policyparameters =list(distancingeffect = sqrt(0.6), hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7) 
#policyparameters =list(distancingeffect = 1, hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7) 
#debug = FALSE
#deltat = 1/10



evolve = function (startstate, startt, stopt, deltat = 1/10, R0 = 2, policyparameters =list(distancingeffect = 0.6, hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7) , debug = FALSE ){ 
    t = startt
    t1 = 1
    nstates = ncol(startstate)
    npops=nrow(startstate)
    PP=startstate[,1]
    nt = stopt - startt
    states = rep(0, npops*nstates*ceiling(nt/deltat))
    states = matrix (states, ceiling(nt/deltat), nstates*npops)
    colnames(states)=as.vector(outer(c('S','E','I','HH','HC','C','R'), c(1:npops), paste, sep="."))
    dim(states)
    head(states)
    times = rep(0, ceiling(nt/deltat))

    currentstate = startstate
    cumstate = currentstate
    cases = rep(0,nrow(states))
    allcases = 0 

    states[t1,] = as.vector(t(currentstate))
    times[t1] = t
    cases[t1] = 0
    
    Qmat=makeQmatrix(policyparameters, startt, stopt,nt, deltat,npops)
    
    
    distancing = 0
    trans = maketransmatrix_MultiPop (t, currentstate, distancing, R0, policyparameters ) 
    
    vals = rep(0, nrow(trans)*ncol(trans))


    while (t<stopt){  
        delta = rep(0, nstates*npops)


        if (debug) { 
            cat ("*******\n")
            cat ("time = " ,t,"\n")
            cat ("current state = ", currentstate,"\n")
        }

        ####distancing = distancingpolicy (t, currentstate, distancing, policyparameters)
        # Time-dependent transition matrix
        # Can change depending on policy
        # Can also modify to have seasonal effects like in Kissler et al.
        distancing=0
        trans = maketransmatrix (t, currentstate, distancing, R0, policyparameters) 

        sd=policyparameters$distancingeffect*Qmat[,ceiling(t*deltat)]
        sd[which(sd==0)]=1
        
        for (popnum in 1:npops){
          
         for (i in 1:nrow(trans)) {
              tmp = trans[i,]

              # For each transition, get the rate. 
              val = tmp[3]

              # Handling the non-linear transition: S->E  = \beta SI
              if ( tmp[1]==1 && tmp[2]==2){ 
                
                  
                
                #val = val * currentstate[,3]*sd*PP
                val = val * currentstate[,3]*sd
                
                  # In the initial few time steps, this allows for some individuals to go from S->E (since the default start state is c(1,0,0,0,0,0,0))    
                #if (t < startt + 3.5 | sum(currentstate[,3])<(0.01/7)){
                if (t < startt + 3.5 | sum(currentstate[popnum,1])<(0.01/7)){
                    val = val + 0.01/7*PP
                  }
                  val = sum(val * currentstate[popnum,1] * sd[popnum])
              } else {
                 val = val * currentstate[popnum,tmp[1]]
              }
              vals[i+(popnum-1)*7] = val        
              delta[tmp[1]+(popnum-1)*7] = delta[tmp[1]+(popnum-1)*7] - val     
              delta[tmp[2]+(popnum-1)*7] = delta[tmp[2]+(popnum-1)*7] + val     
  
              if (tmp[2]==4 || tmp[2]==5){
                  allcases = allcases + val
              }

          }   # trans loop     
        } #end popnum loop
       
        
          if (t1 > nrow(states)){
              states = rbind (states,as.vector(t(currentstate)))
              times = c(times,t)
              cases = c(cases, allcases)
          } else {
              states[t1,] = as.vector(t(currentstate))
              times[t1] = t
              cases[t1] = allcases
          }        
          delta = delta * deltat;
          t = t + deltat
          t1 = t1 + 1
          if (debug) {
            cat ("trans = ", vals,"\n")
           }

          currentstate = currentstate + matrix(delta, nrow=nrow(currentstate),ncol=ncol(currentstate),byrow = T)
          cumstate = cumstate + currentstate
          
        
    } #end while t < stopt
    rownames (states) = times
    return (list(states=states, startt = startt, stopt = stopt, startstate = startstate, trans = trans, times = times, cases = cases)) 
}

