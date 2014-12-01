data(LSAT7)
LSAT7 = as.matrix(LSAT7)
datos = expand.table(LSAT7)
#datos = read.table("file:///home/mirt/ValidaciónTodoElAlgoritmo/Datasets/Test_10_6_1000.csv",sep=" ",header=T)
#datos = read.table("file:///home/mirt/ValidaciónTodoElAlgoritmo/Datasets/Test_10_2_1000.csv",sep=" ",header=T)
#datos = read.table("file:///home/mirt/ValidaciónTodoElAlgoritmo/Datasets/Test_10_10_2000.csv",sep=" ",header=T)
#datos = read.table("file:///home/mirt/ValidaciónTodoElAlgoritmo/Datasets/Test_100_2_10000.csv",sep=" ",header=T)

datos = read.table("file:///home/mirt/Validaciones_Modelos_Principales/3PL/Datasets/Test_10_1_1000.csv",sep=" ",header=T)


inicio = Sys.time()
est = estimacion.Newton(datos)
Sys.time() - inicio
#sink()

#Ajuste Mirt
inicioM = Sys.time()
fit = mirt(data=datos,model=1,itemtype="3PL",verbose = T)
Sys.time() - inicioM
coef.mirt  = unlist(coef(fit))
coef.mirt = coef.mirt[1:(length(coef.mirt)-2)]
coef.mirt = matrix(coef.mirt,ncol = 4,byrow=T)
coef.mirt = coef.mirt[,1:3]

est$zita
coef.mirt

coef.mirt - est$zita

mean(abs(coef.mirt - est$zita))
max(abs(coef.mirt - est$zita))


est$contadorNear
est$ciclos

