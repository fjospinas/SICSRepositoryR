

plot.item(est = est,item = 3,numboot = 1000,alpha = 0.05)


est$zita
parsMirt = list()
fit = mirt(data=datos,model=1,itemtype="3PL",verbose = T,SE = T)
itemplot(fit,1,CE=T)
