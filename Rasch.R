#Valores iniciales Andrade

start.andrade = function(datos){
  I = ncol(datos)
  P = nrow(datos)
  m = 5
  
  #scores
  T = apply(datos,1,sum)
  
  #correlacion biserial
  corr.bis = rep(NA,I) 
  for(i in 1:I){
    corr.bis[i]  = cor(datos[,i],T,method="pearson")    
  } 
  #a inicial
  a.ini = sqrt(corr.bis^2/(1-(corr.bis^2)))
  #proporcion de aciertos
  pi = as.vector(apply(datos,2,sum) / P)
  #a inicial
  b.ini = -(qnorm(pi) / corr.bis)
  #b inicial
  c.ini = 1 / rep(m,I) 
  ini.andrade = matrix(c(a.ini,b.ini,c.ini),ncol=3)
  colnames(ini.andrade) = c("a","b","c")
  ini.andrade
}


patrones = function(datos){
  items = ncol(datos)
  comprim = apply(datos,MARGIN=1,FUN=paste,collapse="/")
  freq = table(comprim)
  pats = names(freq)
  freq = as.vector(freq)
  pats = as.numeric(unlist(lapply(pats,FUN=strsplit,split="/")))
  pats = matrix(pats,ncol=items,byrow=T)
  pats = cbind(pats,freq)
  pats
}

library(mirt)
datos = expand.table(LSAT7)


pats = patrones(datos)
patsSinCol = pats[,-ncol(pats)]
ini.and = start.andrade(datos)

xpuntoi = colSums(patsSinCol)
rv = rowSums(patsSinCol)
xvpunto = rv[2]
b = ini.and[,2]
rv
xpuntoi

ceroTheta = function(theta){
  sum(exp(theta - b) / (1-exp(theta - b))) - xvpunto
}
ceroTheta(2.5)

x = seq(from = -5,to = 5,by = .1)
y = rep(NA,length.out = length(x))

for(i in 1: length(x)){
  y[i] = ceroTheta(x[i],rv[2],ini.and[,2])
}
plot(x,y,type = "l")
uniroot(f = ceroTheta,interval = c(-5,5),rv[2],b = ini.and[,2])

secante = function(f,x1,x2,eps = 10^(-5),num = 1000,...){
  while (abs(x1 - x2) > eps) {
   print("Abs")
    print(abs(x1-x2))
    print("f")
    print(f(x2))
   print("x2")
   print(x2)
    c = x2 - f(x2) * (x2 - x1)/(f(x2) - f(x1))
    x1 = x2
    x2 = c
  }
  x2
}

sal =secante(f = ceroTheta,x1 = -3,x2 = -2,xvpunto = rv[2],b = ini.and[,2])
sal
rm(xvpunto
   )
xvp
