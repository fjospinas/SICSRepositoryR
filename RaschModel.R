library(mirt)

data = read.table(file = "file:///home/mirt/Validaciones_Modelos_Principales/Bloque 1/3PL/Datasets/Test_50_4_5000.csv" ,header = T,sep = " ")

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


func = function(vecB,data){
  score = rowSums(data)
  itemScore = colSums(data)
  scoreUniq = sort(unique(score))
  
  patsSinFrec = cbind(patsSinFrec,rowSums(patsSinFrec))
  agrup = list()
  frecScore = rep(0,length(scoreUniq))
  
  for(i in 1:length(scoreUniq)){
    agrup = append(agrup,list(which(patsSinFrec[,ncol(patsSinFrec)] == scoreUniq[i])))
    frecScore[i] = length(which(rowSums(data) == scoreUniq[i]))
  }
  
  print(frecScore)
  print("--------------")
  print(agrup)
  
  
  gammaR = rep(0,length(agrup))
  for(i in 1:length(agrup)){
    inds = agrup[[i]]
    n = length(inds)
    print(i)
    print(n)
    if(n == 1 ){
      subMat = pats[inds,][-ncol(pats)]
    }else{
      subMat = pats[inds,][,-ncol(pats)]
    }    
    matB = matrix(rep(vecB,n),nrow=n,byrow = T)
    prod = (matB * subMat)
    expo = exp(rowSums(prod))
    gammaR[i] = sum(expo)
  }
  
  print(gammaR)
  
  denom = prod(gammaR ^ frecScore)
  print(gammaR ^ frecScore)
  print(denom)
  
  num =  exp(-sum(itemScore * vecB))
  num /denom
}

func(vecB = b.ini,data = data)


