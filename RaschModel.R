library(mirt)
data = expand.table(LSAT7)

pats = patrones(data)
pats

and = start.andrade(data)
b.ini = and[,2]
b.ini

prob = function(theta, beta){
  1 / (1 + exp(-(theta-beta)))
}

start.theta = function(datos){
  score = rowSums(datos)
  (score - mean(score)) / sd(score) 
}

theta.ini = start.theta(data)



func.max.theta = function(theta,b){
  n = length(b)
  sum = 0
  for(i in 1:n){
    p = prob(theta = theta,beta = b[i])
    sum = sum + (p * (1-p))
  }
  sum
}

func.max.theta(theta.ini[1],b.ini)
salida = nlm(f = func.max.theta,p = theta.ini[1000],b = b.ini)
salida$estimate

length(theta.ini)
