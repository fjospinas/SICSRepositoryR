#Probabilidad
gg = function(a,d, cp,  theta){
  exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-D*(a*theta+ d)))^(-1)
}

personfit = function(est,data){
  scores = scoresEAP(est)  
  nitems = ncol(datos)
  nscores = nrow(scores)
  ninds = nrow(datos)
  
  
  scoresTot = rep(0,nscores)
  for(i in 1:ninds){
    for(j in 1:nscores){
      if(sum(datos[i,] == scores[j,][1:nitems]) == nitems){
        scoresTot[i] = scores[j,][nitems + 1]
      }
    }
  }
  print(scoresTot)
  LL = matrix(0,nrow = ninds,ncol = nitems)
  for(j in 1:ninds){
    for(i in 1:nitems){
      LL[j,i] = gg(a = est$zita[i,1],d = - est$zita[i,2] * est$zita[i,1],cp = qlogis(est$zita[i,3]),theta = scoresTot[j])
    }
  }
  print(LL)
}

personfit(est = est,data = datos)
