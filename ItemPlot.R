library(mvtnorm)

gg = function(a,d, cp,  theta){
  exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-D*(a*theta+ d)))^(-1)
}

est$zita
plot.item = function(zita = est$zita, item,numboot,alpha ){
  nitems = nrow(est$zita)
  x = seq(from = -6,to = 6,by = .1)
  y = sapply(X = x,FUN = gg,a = est$zita[item,1],d = est$zita[item,2],cp = qlogis(est$zita[item,3]))  
  inf = diag(c(est$hess[item,item],est$hess[nitems + item,nitems + item],est$hess[2 * nitems + item,2 * nitems + item]))
  inf[1,2] = inf[2,1] = est$hess[item,nitems + item]
  inf[1,3] = inf[3,1] = est$hess[item,2 * nitems + item]
  inf[2,3] = inf[3,2] = est$hess[nitems + item,2 * nitems + item]
  inf = solve(inf) 
  mean = est$zita[item,]
  boot = rmvnorm(numboot,mean = mean,sigma = inf)
  boot[,3] = ifelse(boot[,3] >= 1,1 - 1e-6,boot[,3])
  boot[,3] = ifelse(boot[,3] <= 0,sqrt(2.2e-16),boot[,3])
  boot[,1] = ifelse(boot[,1] <= 0,sqrt(2.2e-16),boot[,1])
  print("+++++++++++++++++++")
  print(apply(X = boot,MARGIN = 2,max))
  print(apply(X = boot,MARGIN = 2,min))
  print("----------------")
  print(y[1])
  #print(boot)
  #dens = dmvnorm(x = boot,mean = mean,sigma = inf)
  #sortdens = sort(dens)
  #lower <- sortdens[floor(length(dens) * alpha / 2)]
  #upper <- sortdens[ceiling(length(dens) * (1 - alpha / 2))]
  #index = sapply(c(lower,upper),FUN = function(x){which(dens == x)})
  #zitalow = boot[index[1],]
  #zitaup = boot[index[2],]
  #ylow = sapply(X = x,FUN = gg,a = zitalow[1],d = zitalow[2],cp = qlogis(zitalow[3]))
  #yup = sapply(X = x,FUN = gg,a = zitaup[1],d = zitaup[2],cp = qlogis(zitaup[3]))
  envelop = matrix(0,nrow = length(x),ncol = 2)
  for(i in 1:length(x)){
    trace = apply(X = boot,
                  MARGIN = 1,
                  FUN = function(puntoBoot){
                    gg(a =puntoBoot[1], d = puntoBoot[2],cp = qlogis(puntoBoot[3]),theta = x[i])
                    })
    trace = sort(trace)
    lower <- trace[floor(length(trace) * alpha / 2)]
    upper <- trace[ceiling(length(trace) * (1 - alpha / 2))]
    envelop[i,] = c(lower,upper)
    
  }
  plot(x,y,xlim=c(-6,6),ylim=c(0,1),type="l")
  lines(x,envelop[,1],col="red")
  lines(x,envelop[,2],col="blue")

}




plot.item(est$zita,1,1000,0.1)


est$zita
parsMirt = list()
fit = mirt(data=datos,model=1,itemtype="3PL",verbose = T,SE = T)
itemplot(fit,1,CE=T)
