#include <Rcpp.h>
#include <iostream>
//#include <omp.h>
#include <cmath>
//#include <iomanip>

//Función valor absoluto ya que la de cpp no sirve


/*Función que trunca Q
 */


SEXP invMat(SEXP x,bool *llamadoNear){
  using namespace Rcpp ;
  BEGIN_RCPP
  NumericMatrix xx(x);
  NumericMatrix retorno(3,3);

  NumericVector d(3);
 
  //Corrige nan en la matriz
  for(int i = 0 ; i < 3 ; i++){
    for(int k = 0 ; k < 3 ; k++){
      if(::isnan(xx(i,k))){
	xx(i,k) = 0.0; 
      }
    }
  }
  
  Function solve ("solve");
  retorno = solve(xx);
  //Elimina variables y retorna
    
  return(retorno);
  END_RCPP
}

/*Función que calcula la probabilidad
 */
double Pr(int u,double a, double d, double c, double theta){
  double p = 0;
  double z = a * theta + d;
  if(abs(z) > 35){
   z = abs(z) /  z * 35; 
  }
  p = exp(c) / (1+exp(c)) + (1-(exp(c) / (1+exp(c)))) * (1/(1 + exp(-1*(z))));
  if(0.0 == p){
    p = sqrt(2.2e-16);
  }else if(1.0 == p){
    p = 1 - 1e-6;
  }
  if(1 == u){
    return(p);
  }else{
    return(1-p);
  }
}

/*Función que calcula la probabilidad de los moños
 */
double Prm(int u,double a, double d, double theta){
  double p = 0;
  double z = a * theta + d;
  if(abs(z) > 35){
   z = abs(z) /  z * 35; 
  }
  p =  (1/(1 + exp(-1*(z))));
  if(0.0 == p){
    p = sqrt(2.2e-16);
  }else if(1.0 == p){
    p = 1 - 1e-6;
  }
  if(1 == u){
    return(p);
  }else{
    return(1-p);
  }
}

RcppExport SEXP pasoM(SEXP zita,SEXP Rmat,SEXP f, SEXP ptcuad, SEXP wcuad, SEXP iniVals){
  using namespace Rcpp ;
  BEGIN_RCPP
  //Conversion de entradas
  NumericMatrix xzita(zita);
  NumericMatrix xzitaAnt = clone(xzita);
  NumericMatrix xiniVals(iniVals);
  NumericMatrix xRmat(Rmat);
  NumericVector xf(f);
  NumericVector xptcuad(ptcuad);
  NumericVector xwcuad(wcuad);
  
  //Constantes
  const int items = xzita.ncol();
  const int numCuads = xptcuad.length();
  
  //Variables utilitarias
  int i, k, j, m;
  //Calculo de g y h
  NumericMatrix g(3,items);
  NumericVector h(3);
  NumericMatrix retorno(3,items);
  bool *llamadoNear = new bool;
  *llamadoNear = 0;
  
  for(i = 0 ; i < items ; i++){
   double *suma = new double[3];
   //Inicializa suma a cero
   for(j = 0 ; j < 3 ; j++){
     suma[j] = 0;
   }
   for(k = 0 ; k < numCuads ; k++){
     double p = Pr(1,xzita(0,i), xzita(1,i) ,xzita(2,i), xptcuad(k));
     double pm = Prm(1,xzita(0,i), xzita(1,i), xptcuad(k));
     double aux = (xRmat(k,i) -  xf[k] * p) * ((pm * (1 - pm)) / (p * (1 - p)));
     h[0] = xptcuad[k]*(1/(1+exp(xzita(2,i))));
     h[1] = (1/(1+exp(xzita(2,i))));
     h[2] = pow((1/(1+exp(xzita(2,i)))),2) * exp(xzita(2,i)) / pm;
     for(j = 0 ; j < 3 ; j++){
      suma[j] = suma[j] + (aux * h[j]); 
     }
   }
   for(j = 0 ; j < 3 ; j++){
     g(j,i) = suma[j]; 
   }
   delete[] suma;
  }
  
  //Actualización de zita
  for(i = 0 ; i < items ; i++){
    NumericMatrix suma(3,3);
    //Inicializa suma
    for(j = 0 ; j < 3 ; j++){
     for(k = 0 ; k < 3 ; k++){
      suma(j,k) = 0.0; 
     }
    }
    for(k = 0 ; k < numCuads ; k++){
      double p = Pr(1,xzita(0,i), xzita(1,i) ,xzita(2,i), xptcuad(k));
      double q = 1.0 - p;
      double pm = Prm(1,xzita(0,i), xzita(1,i), xptcuad(k));
      double qm = 1.0 - pm;
      
      //Gradiente probabilidad
      h = NumericVector(3);
      h[0] = xptcuad[k] * (1 / (1+exp(xzita(2,i)))) * pm * qm;
      h[1] = (1/(1+exp(xzita(2,i)))) * pm * qm;
      h[2] = pow((1/(1+exp(xzita(2,i)))),2) * exp(xzita(2,i)) * qm;
      
      NumericMatrix hess(3,3);
      hess(0,0) = (pow(xptcuad[k],2))/(1+exp(xzita(2,i)))*(1-2*pm);
      hess(0,1) = hess(1,0) = ((xptcuad[k])/(1+exp(xzita(2,i)))*(1-2*pm));
      hess(0,2) = hess(2,0) = ((-exp(xzita(2,i)) * xptcuad[k])/pow((1+exp(xzita(2,i))),2));
      hess(1,1) = ((1-2*pm)/(1+exp(xzita(2,i))));
      hess(1,2) = hess(2,1) = -(exp(xzita(2,i)) / pow((1+exp(xzita(2,i))),2));
      hess(2,2) = ((exp(xzita(2,i)) * (1-exp(xzita(2,i))))/(pm * pow((1+exp(xzita(2,i))),3)));
      for(j = 0 ; j < 3 ; j++){
	for(int m = 0 ; m < 3 ; m++){
	  hess(j,m) = hess(j,m) * pm * qm;
	}
      }
            
      //Calcula la suma final en multiples pasos
      double aux1 = (xRmat(k,i) - xf[k] * p) * (1 / (p * q));
      double aux2 = (xRmat(k,i) - 2*p*xRmat(k,i) + xf[k] * pow(p,2)) * (pow((1 / (p * q)),2));
      
      //producto h %*% t(h)
      NumericMatrix prodh(3,3);
      for(m = 0 ; m < 3 ; m++){
	for(j = 0 ; j < 3 ; j++){
	  prodh(m,j) = h(m) * h(j);
	}
      }
      
      for(m = 0 ; m < 3 ; m++){
	for(j = 0 ; j < 3 ; j++){
	  suma(m,j) = suma(m,j) + aux1 * hess(m,j) - aux2 * prodh(m,j);
	}
      }
    }      
    
    // zita[,i] = zita.ant[,i] - inv.mat(delta.zita[,,i]) %*% g[,i]
    NumericMatrix invSuma(3,3);
    //Function solve ("solve");
    //invSuma = solve(suma);
    invSuma = invMat(suma,llamadoNear);
    NumericVector prodInvG(3);
    for(j = 0 ; j < 3 ; j++){
      double sum = 0;
     for(m = 0 ; m < 3 ; m++){
       sum = sum + invSuma(j,m) * g(m,i);        
     }
     prodInvG(j) = sum;
    }
    
    for(j = 0 ; j < 3 ; j++){
      retorno(j,i) = xzitaAnt(j,i) - prodInvG(j);
    }
    
  }
  
  //Cotas para los valores
  for(i = 0 ; i < items ; i++){
    if(std::abs(retorno(0,i)) > 10) retorno(0,i) = xiniVals(0,i);
    if(std::abs(retorno(1,i)) > 40) retorno(1,i) = xiniVals(1,i);
    if(std::abs(retorno(2,i)) > 600) retorno(2,i) = xiniVals(2,i);
  }
  //Retorna zita
  IntegerVector llamado(1);
  if(*llamadoNear){
    llamado(0) = 1;
  }else{
    llamado(0) = 0;
  }
  
  
  delete llamadoNear;
  return  List::create(_["zita"] = retorno,_["llamadoNear"] = llamado);
  END_RCPP
}

//Función de logverosimilitud
RcppExport SEXP Loglik(SEXP zita,SEXP Rmat,SEXP f, SEXP ptcuad){
  using namespace Rcpp ;
  BEGIN_RCPP
  NumericVector xzita(zita);
  NumericMatrix xRmat(Rmat);
  NumericVector xf(f);
  NumericVector xptcuad(ptcuad);
  
  const int items = xzita.length() / 3;
  const int numCuads = xptcuad.length();
  
  double suma = 0.0;
  for(int k = 0 ; k < numCuads ; k++){
    for(int i = 0 ; i < items ; i++){
      double pki = Pr(1, xzita[i], xzita[items + i], xzita[2 * items + i], xptcuad[k]);
      double qki = 1.0 - pki;
      suma = suma + (xRmat(k,i) * log(pki) + (xf[k] - xRmat(k,i)) * log(qki));
    }    
  }
  NumericVector ret(1);
  ret(0) = -suma;
  return ret;
  END_RCPP
}

//Gradiente logverosimilitud
RcppExport SEXP grad(SEXP zita,SEXP Rmat,SEXP f, SEXP ptcuad){
  using namespace Rcpp ;
  BEGIN_RCPP
  NumericVector xzita(zita);
  NumericMatrix xRmat(Rmat);
  NumericVector xf(f);
  NumericVector xptcuad(ptcuad);
  
  const int items = xzita.length() / 3;
  const int numCuads = xptcuad.length();
  
  NumericMatrix g(3,items);
  NumericVector h(3);
  NumericMatrix retorno(3,items);
  bool *llamadoNear = new bool;
  *llamadoNear = 0;
  
  for(int i = 0 ; i < items ; i++){
   double *suma = new double[3];
   //Inicializa suma a cero
   for(int j = 0 ; j < 3 ; j++){
     suma[j] = 0;
   }
   for(int k = 0 ; k < numCuads ; k++){
     double p = Pr(1,xzita[i], xzita[items + i] ,xzita[2 * items + i], xptcuad(k));
     double pm = Prm(1,xzita[i] , xzita[items + i], xptcuad(k));
     double aux = (xRmat(k,i) -  xf[k] * p) * ((pm * (1 - pm)) / (p * (1 - p)));
     h[0] = xptcuad[k]*(1/(1+exp(xzita[2 * items + i])));
     h[1] = (1/(1+exp(xzita[2 * items + i])));
     h[2] = pow((1/(1+exp(xzita[2 * items + i]))),2) * exp(xzita[2 * items + i]) / pm;
     for(int j = 0 ; j < 3 ; j++){
      suma[j] = suma[j] + (aux * h[j]); 
     }
   }
   for(int j = 0 ; j < 3 ; j++){
     g(j,i) = suma[j]; 
   }
   delete[] suma;
  }
  
  NumericVector grad(items * 3);
  for(int j = 0 ; j < 3 ; j++){
   for(int i = 0 ; i < items ; i++){
    grad(j*items + i) = -g(j,i); 
   }
  }
  return grad;
  END_RCPP
}
