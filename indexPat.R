#Script que indexa los datos a su respectivo patron.
#es necesario para retornar las habilidades expandidas

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

pats = patrones(data)

indexPat = function(data,pats){
  comprimData = apply(data,MARGIN=1,FUN=paste,collapse="/")
  comprimPats = apply(pats[,-ncol(pats)],MARGIN=1,FUN=paste,collapse="/")
  index = lapply(comprimPats,FUN = function(x) {which(x == comprimData)})
  index
}

index = indexPat(data,pats)
index
pats
