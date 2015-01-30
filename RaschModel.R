library(mirt)
data = read.table(file = "file:///home/mirt/Validaciones_Modelos_Principales/Bloque 1/3PL/Datasets/Test_10_1_1000.csv" ,header = T,sep = " ")

pats = patrones(data)
patsSinFrec = pats[,-ncol(pats)]

and = start.andrade(data)
b.ini = and[,2]
b.ini

start.theta = function(datos){
  score = rowSums(datos)
  (score - mean(score)) / sd(score) 
}

theta.ini = start.theta(data)
theta.ini

score = rowSums(data)
score

scoreUniq = sort(unique(score))
scoreUniq

patsSinFrec = cbind(patsSinFrec,rowSums(patsSinFrec))


agrup = list()
frecScore = rep(0,length(scoreUniq))

for(i in 1:length(scoreUniq)){
  agrup = append(agrup,list(which(patsSinFrec[,ncol(patsSinFrec)] == scoreUniq[i])))
  frecScore[i] = length(which(rowSums(data) == scoreUniq[i]))
}

gammaR = rep(0,length(agrup))
for(i in 1:length(agrup)){
  inds = agrup[[i]]
  n = length(inds)
  subMat = pats[inds,][,-ncol(pats)]
  matB = matrix(rep(b.ini,n),nrow=n,byrow = T)
  prod = (matB * subMat)
  expo = exp(rowSums(prod))
  gammaR[i] = sum(expo)
}

denom = prod(gammaR ^ frecScore)

