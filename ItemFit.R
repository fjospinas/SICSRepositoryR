#Item.fit de SICS
itemFit = function(est,G = 10,FUN = median,p.val.sim = FALSE,boot.num = 100){
  pats = est$pats
  zita = est$zita
  scores = scoresEAP(est)
  nitems = nrow(zita)
  theta = scores[,nitems + 1 ]  

  
  frec = pats[,ncol(pats)]
  patsSinFrec = pats[,-ncol(pats)]
  groups = quantile(theta,seq(0, 1, length.out = G + 1))
  print("groups")
  print(groups)
  groups[1] = groups[1] - 0.1
  print("groups")
  print(groups)
  groups[G + 1] = groups[G + 1] + 0.1
  print("groups")
  print(groups)
  groups.ind = findInterval(theta,groups)
  print("gi")
  print(groups.ind)
  print("ugi")
  print(sort(unique(groups.ind)))
  groups.ind = factor(groups.ind, levels = sort(unique(groups.ind)))
  print("gi")
  print(groups.ind)
  print("gi")
  print(groups.ind)
  thetaG = tapply(rep(theta, frec), rep(groups.ind, frec), FUN = FUN)
  print("groups")
  print(groups)
  print("gi")
  print(groups.ind)
  print("tg")
  print(thetaG)
  pr = matrix(NA,ncol = nitems,nrow = G)
  for(i in 1:G){
    for(j in 1:nitems){
      pr[i,j] = gg(a = zita[j,1],d = zita[j,2],cp = qlogis(zita[j,3]),theta = thetaG[i])
    }
  }
  Nj = as.vector(tapply(frec, groups.ind, sum))
  print("nj")
  print(Nj)
  
  Obs = rowsum(frec * patsSinFrec, groups.ind, reorder = FALSE)
  Obs2 = rowsum(frec * patsSinFrec, groups.ind, reorder = TRUE)
  print("obs2")
  print(Obs2)
  
  print("pr")
  print(pr)
  
  X2 = numeric(nitems)
  for(i in 1:nitems){
    for(j in 1:G){
      P = pr[j,i]
      P = c(1-P,P)
      r = c(Nj[j] - Obs2[j,i],Obs2[j,i])
      print("------------")
      print(paste(j,",",i))
      print(P)
      print(r)
      X2[i] = X2[i] + sum((r - Nj[j]*P)^2 / Nj[j]*P)
    }
  }
  print("X2")
  print(X2)
  #chi.square = Nj * (Obs - pr)^2/(pr * (1 - pr))
  #Tobs = colSums(chi.square, na.rm = TRUE)
  #df = G - 3 
  #pvals <- pchisq(Tobs, df = df, lower.tail = FALSE)
  #salida = matrix(c(Tobs,pvals),ncol= 2)
  salida
}

a=itemFit(est,G=10)
