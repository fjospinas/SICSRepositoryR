#include <Rcpp.h>
#include <iostream>
#include <omp.h>

//Función que calcula la probabilidad
double Pr(int u,double a, double b, double c, double theta){
  double p = 0;
  p = exp(c) / (1+exp(c)) + (1-(exp(c) / (1+exp(c)))) * (1/(1 + exp(-1*(a* theta + b))));
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

//Fución que calcula R y f y los retorna como una lista a R
//Se invoca desde R
RcppExport SEXP calculoRF2(SEXP zita, SEXP theta,SEXP A,SEXP pats){
  using namespace Rcpp ;
  BEGIN_RCPP
  const NumericMatrix xzita(zita);
  const NumericVector xtheta(theta);  
  const NumericVector xA(A);
  const NumericMatrix xpats(pats);
  
  const int items = xzita.ncol();
  const int numCuads = xtheta.length();
  const int nPats = xpats.nrow();  
  
  NumericMatrix Rmat(numCuads,items);
  NumericVector Fmat(numCuads);
  
  double *faux = new double[numCuads];
  double sum = 0.0;
  
  int k,i,j;
  
  //Calcula R y F
  for(j = 0; j < nPats ; j++){
    //Inicializa faux
    for(k = 0; k < numCuads;k++ ){
      faux[k] = 1; 
    } //fin for    
    
    //calcula g*() para todos los k
    for(k = 0; k < numCuads; k++){
      for(i = 0; i < items; i++){
	faux[k] = faux[k] * Pr((int) xpats(j,i),xzita(0,i), xzita(1,i) ,xzita(2,i) ,xtheta[k]);
      } //fin for
      faux[k] = faux[k] * xA[k];
    } // fin for
    
    //Calcula el denominador de g*()
    sum = 0.0;
    for(k = 0; k < numCuads; k++){
     sum += faux[k]; 
    } // fin for
    
    for(k = 0; k < numCuads; k++){
     faux[k] = faux[k] / sum;
     
     //Se multiplica por la frecuencia del patrón y se almacena en Fmat para retornar
     faux[k] = faux[k] * xpats(j,items);
     Fmat[k] = Fmat[k] + faux[k];
     
     //Se calcula Rmat y se retorna
     for(i = 0 ; i < items ; i++){
       if((int) xpats(j,i)){
	Rmat(k,i) = Rmat(k,i) + faux[k]; 
       } // fin if
     } // fin for
    } // fin for
  } // fin for
  
  //Se eliminan las vaiables y se retornan R y fvec
  delete[] faux;
  return List::create(_["R"] = Rmat,_["fvec"] = Fmat);
  END_RCPP
} 