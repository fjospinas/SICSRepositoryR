#Probabilidad
gg = function(a,d, cp,  theta){
  exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-D*(a*theta+ d)))^(-1)
}

itemFit2 = function(est,data){
  scores = scoresEAP(est)  
  nitems = ncol(data)
  nscores = nrow(scores)
  ninds = nrow(data)
  
  
  scoresTot = rep(0,nscores)
  for(i in 1:ninds){
    for(j in 1:nscores){
      if(sum(data[i,] == scores[j,][1:nitems]) == nitems){
        scoresTot[i] = scores[j,][nitems + 1]
      }
    }
  }
  
  LL = P = Q = matrix(0,nrow = ninds,ncol = nitems)
  for(j in 1:ninds){
    for(i in 1:nitems){
      P[j,i] = gg(a = est$zita[i,1],d = est$zita[i,2],cp = qlogis(est$zita[i,3]),theta = scoresTot[j])
    }
  }
  
  print(dim(LL))
  Q = 1 - P
  LL[data == 1] = P[data == 1]
  LL[data == 0] = Q[data == 0]
  LL = colSums(log(LL))
  print(LL)
  
  mu = sigmaCuad = rep(0,nitems)  
  for( i in 1:nitems){
    Pi = cbind(P[,i],Q[,i])
    logPi = log(Pi)
    mu[i] = sum(Pi * logPi)    
    #sigmaCuad = sigmaCuad + Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2)
    sigmaCuad[i] = sum(Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
    
  }
  print("mu")
  print(mu)
  print(sigmaCuad)
  Z3 = (LL - mu) / sqrt(sigmaCuad)
  Z3
}


z=itemFit2(est = est,data = data)
z
