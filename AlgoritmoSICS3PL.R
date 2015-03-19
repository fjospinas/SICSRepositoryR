rm(list = ls(all = TRUE))
library(mirt)
library(Matrix)
library(numDeriv)
library(optimx)

#-fopenmp

setwd("/home/mirt/Git/GrupoSICS/dev/SICSRepositoryR/")
system("PKG_CPPFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'` PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` R CMD SHLIB  pasoe3.cpp")
system("PKG_CPPFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'` PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` R CMD SHLIB  pasom3.cpp")
dyn.load("pasoe3.so")
dyn.load("pasom3.so")

###########
# PRUEBAS #
###########

D = 1

##################################################
# 3.Funcion para inicializar valores con Andrade #
##################################################

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

#función de puntos de cuadratura y pesos
library(statmod)
Cuad = gauss.quad(n=41,"hermite")
pt.cuad = Cuad[[1]] * sqrt(2)
pt.cuad = - (pt.cuad)
#pt.cuad
w.cuad = Cuad[[2]]  /  sqrt(pi)

pt.cuad = read.table("/home/mirt/Trabajo IRT/Algoritmo SICS/PWcuad.csv",dec=".",sep = " ",header = T)
w.cuad = pt.cuad[,2]
pt.cuad = pt.cuad[,1]


#Valores Iniciales Andrade
#and = t(start.andrade(datos))
#and[2,] = -and[2,] * and[1,]
#and[3,] = qlogis(and[3,])

#valores iniciales MIRT
inicio.mirt = function(datos){
  ini.mirt = mirt(data=datos,model=1,itemtype="3PL",pars="values")
  ini.mirt = ini.mirt$value
  ini.mirt = ini.mirt[1:(length(ini.mirt) -2)]
  ini.mirt = t(matrix(ini.mirt,ncol=4,byrow=T)[,c(1,2,3)])
  ini.mirt[3,] = qlogis(ini.mirt[3,])
  and = ini.mirt
  and
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

#Probabilidad
gg = function(a,d, cp,  theta){
  exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-D*(a*theta+ d)))^(-1)
}


#Log verosimilitud
LL = function(zita.vec,R,fvec,pt.cuad,nitems){
  suma = 0
  for(k in 1:41){
    for(i in 1:nitems){
      rki = R[k,i]
      fki = fvec[k]
      a = zita.vec[i]
      d = zita.vec[nitems + i]
      c = zita.vec[2*nitems + i]
      pki = gg(a=a,d=d,cp=c,theta=pt.cuad[k])
      qki = 1 - pki 
      suma = suma + (rki*log(pki)+(fki-rki)*
                       log(qki))
      
    }
  }
  -suma
}

LL2 = function(zita.vec,R,fvec,pt.cuad,nitems,and){
   .Call("Loglik",zita.vec,R,fvec,pt.cuad)
}

gradLoglik = function(zita.vec,R,fvec,pt.cuad,nitems,and){
  .Call("grad",zita.vec,R,fvec,pt.cuad)
}


inicio = Sys.time()
###################
# Algoritmos SICS #
###################
estimacion.Newton = function(datos){
  and = inicio.mirt(datos)
  and.copia = and
  pats = patrones(datos)
  npats = nrow(pats)
  nitems = ncol(pats) - 1
  listaOptimX = list()
  
  zita.ant = zita = and
  seguir = TRUE
  mm = 0
  contadorNear = 0
  while(seguir){
    inicio.ciclo = Sys.time()
    mm = mm +1
  
    ##########
    # Paso E #
    ##########
    
    #inicioE = Sys.time()  
    RyF=.Call("calculoRF2",and,pt.cuad,w.cuad,pats)
    #te = Sys.time() - inicioE
    #print("*********Tiempo E")
    #print(te)
    R = RyF$R
    fvec = RyF$fvec
  
    ##########
    # Paso M #
    ##########
    
    print("Entra a optim")
    zita.vec = as.vector(t(zita))
    #opt = optim(par=zita.vec,fn=LL,method="BFGS",R=R,fvec=fvec,pt.cuad=pt.cuad,nitems = nitems,control=list(maxit=10))
    opt = optim(par=zita.vec,fn=LL2,gr=gradLoglik,method= "BFGS",R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems,and=and,control=list(maxit=20),hessian = T)
    #optx = optimx(par=zita.vec,fn=LL2,gr=gradLoglik,itnmax = 20,control=list(all.methods=TRUE, save.failures=TRUE, trace=0),R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems,and=and)
    optx = optimx(par=zita.vec,fn=LL2,gr=gradLoglik,itnmax = 20,control=list(all.methods=TRUE, save.failures=TRUE, trace=0),R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems,and=and)
    listaOptimX = append(listaOptimX,list(optx))
    #opt = optim(par=zita.vec,fn=LL2,method= "L-BFGS-B",R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems,and=and,control=list(maxit=10))
                #,lower = c(rep(-10,10),rep(-40,10),rep(-600,10)),upper = c(rep(10,10),rep(40,10),rep(600,10)))
    #opt = vmmin(fr=LL,x=zita.vec,R=R,fvec=fvec,pt.cuad=pt.cuad,nitems = nitems)
    contadorNear = contadorNear + 1
    zita = matrix(opt$par,ncol=nitems,byrow=T)
    hess = opt$hessian
    
    zita[1,] = ifelse(abs(zita[1,]) > 10, and[1,], zita[1,])
    zita[2,] = ifelse(abs(zita[2,]) > 40, and[2,], zita[2,])
    zita[3,] = ifelse(abs(zita[3,]) > 600, and[3,], zita[3,])
    ##Imprime salidas
    
    #print(mean(abs((zita - zita.ant)/zita.ant)))
    print(paste("Fin ciclo: ", mm, " Convergencia: ", max(abs((zita - zita.ant)))," Tiempo Ciclo: ",Sys.time() - inicio.ciclo))
    if(max(abs((zita - zita.ant))) < 10^(-3)){seguir = FALSE}
    and = zita.ant = zita
    if(mm > 200){
      print(paste("El algoritmo superó los ",mm - 1," ciclos",sep=""))
      break()
      #stop(paste("El algoritmo superó los ",mm - 1," ciclos",sep=""))
    }
  } #fin while
  zita[3,] = plogis(zita[3,])
  zita = t(zita)
  list(zita=zita,contadorNear=contadorNear,ciclos = mm,pats = pats,hess = hess,listaOptimX = listaOptimX)
}
gcc = NULL
#sink("/home/mirt/Trabajo IRT/Algoritmo SICS/SalidaAUX.txt")
