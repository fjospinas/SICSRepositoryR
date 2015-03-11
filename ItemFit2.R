#Probabilidad
M3PL <- function(zita,theta){
  exp(zita[3])/(1+exp(zita[3]))+ (1-(exp(zita[3])/(1+exp(zita[3]))))*(1 + exp(-D*(zita[1]*theta+ zita[2])))^(-1)
}

#Función que indexa patrones para expandirlos
indexPat = function(data,pats){
  comprimData = apply(data,MARGIN=1,FUN=paste,collapse="/")
  comprimPats = apply(pats[,1:ncol(data)],MARGIN=1,FUN=paste,collapse="/")
  index = lapply(comprimPats,FUN = function(x) {which(x == comprimData)})
  index
}

#función  que calcula item fit basado en Z3 
itemFit2 = function(est,data){
  zita  = est$zita
  zita[,3] = qlogis(zita[,3])
  scores = scoresEAP(est)  
  nitems = ncol(data)
  nscores = nrow(scores)
  ninds = nrow(data)
  
  #Expansión de patrones sobre los datos originales
  index = indexPat(data,est$pats)  
  scoresTot = numeric(nrow(data))
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] = scores[mm,ncol(data) +1]        
  }
  
  #Matriz de probabilidad
  P = apply(zita,1,M3PL,scoresTot)
  
  #Calculo de logverosimilitud
  LL = matrix(0,ncol = ncol(P),nrow = nrow(P))
  LL[data == 1] = P[data == 1]
  LL[data == 0] = 1 - P[data == 0]
  LL = colSums(log(LL))
  
  #Calculo de estimado Z3
  mu = sigmaCuad = rep(0,nitems)  
  for( i in 1:nitems){
    Pi = cbind(P[,i],1 - P[,i])
    logPi = log(Pi)
    mu[i] = sum(Pi * logPi)    
    #sigmaCuad = sigmaCuad + Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2)
    sigmaCuad[i] = sum(Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
    
  }
  Z3 = (LL - mu) / sqrt(sigmaCuad)
  Z3
}

z=itemFit2(est = est,data = data)