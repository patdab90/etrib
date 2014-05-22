etriutils.addVariables <- function(constr, names){
  ## 
  w <- matrix(0, ncol=length(names), nrow=nrow(constr))
  colnames(w) <- names
  constr <- cbind(constr, w)
  return(constr)
}

etriutils.getVariablesCount <- function(constr){
  ## 
  return(ncol(constr))
}