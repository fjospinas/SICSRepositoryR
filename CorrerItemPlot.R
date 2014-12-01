

plot.item(est,1,1000,0.05)


est$zita
parsMirt = list()
fit = mirt(data=datos,model=1,itemtype="3PL",verbose = T,SE = T)
itemplot(fit,1,CE=T)
