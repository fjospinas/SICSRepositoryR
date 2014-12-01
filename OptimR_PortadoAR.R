##Función de optimización cuasi-newton
vmmin = function(fr,x,R,fvec,pt.cuad,nitems){
  #Constantes
  stepredn = 0.2
  acctol = 0.0001
  reltest = 10
  reltol = 1.490116e-08
  n = length(x)
  
  #variables
  accpoint = FALSE
  B = matrix(0.0,n,n)
  c = g = rep(0.0,n)
  count = funcount = gradcount = 0
  D1 = D2 = 0.0
  f = 0.0
  gradproj = 0.0
  ilast = 0
  notcomp = FALSE
  s = 0.0  
  steplength = 1
  Fmin = 0.0
  X = rep(0.0,n)
  Bvec = x
  t = rep(0.0,n)
  
  #Entorno para calcular el gradiente
  myenv <- new.env()
  
  
  #Acá empieza el programa a calcular
  fail = FALSE
  f = fr(x,R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems)
  notcomp = is.infinite(f)
  if(notcomp){
    print(x)
    stop("La funcion no es evaluable en los valores iniciales")
  }else{
    #proceso de minimización
    Fmin = f
    funcount = gradcount = 1
    assign("x", x, envir = myenv)
    assign("R", R, envir = myenv)
    assign("fvec", fvec, envir = myenv)
    assign("pt.cuad", pt.cuad, envir = myenv)
    g = numericDeriv(quote(fr(x,R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems)), c("x"), myenv)
    g = as.vector(attr(g,"gradient"))
    ilast = gradcount
    endwhile = TRUE
    while(endwhile){
      if(ilast == gradcount){ #inicializa B
        B = diag(1,n)
      } #fin if(ilast == gradcount)
      X = Bvec
      c = g
      g.mat = matrix(rep(g,n),ncol=n,byrow=T)
      t = -rowSums((B*g.mat))
      gradproj = sum(t*g)
      if(gradproj < 0){ #se realiza la busqueda lineal
        steplength = 1
        accpoint = FALSE
        seguir = TRUE
        while(seguir){
          count = 0
          Bvec = X + steplength * t
          count = sum((reltest + X) == (reltest + Bvec))
          if(count < n){ #test de convergencia principal
            f = fr(Bvec,R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems)
            notcomp = is.infinite(f) || is.nan(f)
            funcount = funcount + 1
            accpoint = (!notcomp) && (f <= Fmin + gradproj * steplength * acctol)
            if(!accpoint){ #si el punto no es aceptable
              steplength = steplength * stepredn
            } #fin if(!accpoint)
          } #fin if(count < n)     
          if((count == n) || accpoint){ seguir = FALSE }
        } #Fin while(!accpoint)
        flag = abs(f - Fmin) > reltol * (abs(Fmin) + reltol)
        if(!flag){
          count = n 
          Fmin = f
        }
        if(count < n){
          Fmin = f
          assign("x", Bvec, envir = myenv)
          g = numericDeriv(quote(fr(x,R=R,fvec=fvec,pt.cuad=pt.cuad,nitems=nitems)), c("x"), myenv)
          g = as.vector(attr(g,"gradient"))
          gradcount = gradcount + 1
          t = steplength * t
          c = g - c
          D1 = sum(t * c)
          if(D1 > 0){ #si la actualización de la dirección es posible
            c.mat = matrix(rep(c,n),ncol=n,byrow=T)
            s = rowSums(B*c.mat)
            X = s
            D2 = sum(s*c)
            D2 = 1L + D2 / D1
            for(i in 1:n){
              for(j in 1:n){
                B[i,j] = B[i,j] - (t[i]*X[j]+X[i]*t[j]-D2*t[i]*t[j])/D1
              } #Fin for(j in 1:n)
            } #Fin for(i in 1:n)
          }else{ #Update no posible
            print("La actualización de la dirección no es posible. Se reinicia B")
            ilast = gradcount
          }
        }else{ # else if(count < n)
          if(ilast < gradcount){
            count = 0
            ilast = gradcount
          } #Fin if(ilast < gradcount)
        }                        
      }else{ #Else if(gradproj < 0)
        ilast = gradcount 
        count = 0
      } #Fin  if(gradproj < 0)
      if(count == n && ilast == gradcount) { endwhile = FALSE }
    } #Fin while(endwhile)
  } #fin if(notcopm)
  list(par = Bvec,
       value = f,
       func.count = funcount,
       grad.count = gradcount)
}
